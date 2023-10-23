'use strict';

import * as net from 'net';
import winreg = require('winreg');
import * as vscode from 'vscode';
import {
    LanguageClient, LanguageClientOptions, StreamInfo, ServerOptions, MessageSignature
} from 'vscode-languageclient/node';
import {
    Disposable, Uri, TextDocument, OutputChannel, workspace, window, Terminal, TerminalOptions, env
} from 'vscode';
import { spawn } from 'child_process';

import { extensionContext } from './extension';
import * as completion from './completion';
import { EventEmitter } from 'events';
import { isDeepStrictEqual } from 'util';
import * as path from 'path';
import { connected } from 'process';
import { Server } from 'http';

let omlTerminal: Terminal | undefined;

export class LanguageConnection implements Disposable {

    private socket: any;
    private omlPath: string;
    private exeName: string;
    private omlRoot: string;
    private omlThirdParty: string;
    private platForm: string;
    private compData: completion.CompletionData;
    private outputChannel: OutputChannel;
    private client: LanguageClient;
    private emitter: EventEmitter;
    private cbmap: Map<string, any>;
    private isRunning: boolean;

    constructor() {
        this.socket = '';
        this.omlPath = '';
        this.exeName = '';
        this.omlRoot = '';
        this.omlThirdParty = '';
        this.platForm = '';
        this.emitter = new EventEmitter();
        this.cbmap = new Map;
        this.isRunning = false;
        this.client = {} as LanguageClient;
        this.outputChannel = {} as OutputChannel;
        this.setupEnv();
        this.startConnection();
        this.compData = new completion.CompletionData('', '');
        vscode.languages.registerCompletionItemProvider(['oml'], new completion.CompletionItemProvider(this.compData), '');
        const signTriggerCharacters = ['(', ','];
        vscode.languages.registerSignatureHelpProvider(['oml'], new completion.SignatuareHelpProvider(), ...signTriggerCharacters);
        extensionContext.subscriptions.push(vscode.window.onDidCloseTerminal(this.deleteTerminal, this));
    }

    public dispose(): Thenable<void> {
        return this.closeConnection();
    }

    private setupEnv() {

        const config = vscode.workspace.getConfiguration('OML');
        let exePath = config.get<string | undefined>('OML_EXE');
        if (exePath !== undefined) {
            this.omlPath = exePath;
        }
        else {
            this.omlPath = String(process.env['OML_EXE']);
        }
        this.omlRoot = path.join(this.omlPath, '../../../../');

        this.platForm = "win64";
        if (process.platform !== 'win32') {
            this.platForm = "linux64";
        }
        let rootDir = this.omlRoot;
        this.exeName = path.basename(this.omlPath);
        if (this.exeName === "Compose.exe") {
            process.env['HW_ROOTDIR'] = rootDir;
            process.env['ALTAIR_HOME'] = rootDir;
            process.env['HW_UNITY_ROOTDIR'] = path.join(rootDir, 'hwx');
            process.env['OML_APPDIR'] = path.join(rootDir, 'hwx');
            process.env['OML_HELP'] = path.join(rootDir, 'hwx/help/compose/help/en_us/topics/reference/oml_language/');
            process.env['HW_FRAMEWORK'] = path.join(rootDir, '../common/framework/', this.platForm, 'hwx');
            process.env['QT_QPA_PLATFORM_PLUGIN_PATH'] = path.join(String(process.env['HW_FRAMEWORK']), 'bin', this.platForm, 'platforms');
        }
        else if (this.exeName === "omlconsole.exe") {
            process.env['OML_HELP'] = path.join(rootDir, 'help/win/en/topics/reference/oml_language/');
            process.env['OML_THIRDPARTY'] = path.join(rootDir, 'third_party');
            process.env['OML_INTEL_COMP'] = config.get('OML_INTEL_COMP');
            process.env['OML_INTEL_MKL'] = config.get('OML_INTEL_MKL');
            process.env['OML_FFTW'] = config.get('OML_FFTW');
            process.env['OML_MATIO'] = config.get('OML_MATIO');
            process.env['OML_HDF'] = config.get('OML_HDF');
            process.env['OML_QHULL'] = config.get('OML_QHULL');
            this.omlThirdParty = String(process.env['OML_THIRDPARTY']);
        }
        else {
            window.showErrorMessage('OML executable not found, pls check OML_EXE setting.');
        }
    }
    public isConnectionAlive() {
        return this.isRunning;
    }

    public addListener(name: string, callback: any) {
        if (this.cbmap.has(name) === true) {
            this.emitter.removeListener(name, this.cbmap.get(name));
        }
        this.cbmap.set(name, callback);
        this.emitter.addListener(name, callback);
    }
    private setClient(client: LanguageClient) {
        this.client = client;
    }
    private async recvData(data: Buffer) {
        let str = (data.toString()).replace(/}{/g, '},{'); // segregate multiple messages
        str = str.replace(/\n/g, '\\n'); // hack for new lines
        let jObjs = JSON.parse('[' + str + ']');
        //let jObjs = JSON.parse('[' + (data.toString()).replace(/}{/g, '},{') + ']');
        for (var i = 0; i < jObjs.length; i++) {
            if ('' === jObjs[i]) {
                continue;
            }
            switch (jObjs[i].type) {
                case 'status': {
                    this.isRunning = jObjs[i].data === 'ready';
                    break;
                }
                case 'keywordslist': {
                    this.compData.setKeywords(jObjs[i].data);
                    break;
                }
                case 'functionslist': {
                    this.compData.setFunctions(jObjs[i].data);
                    break;
                }
                case 'signatures': {
                    const signs = jObjs[i].data;
                    this.emitter.emit('signature', signs);
                }
            }
        }
    }

    public async sendRequest(request: string) {
        if (!this.isConnectionAlive()) {
            await this.restartConnection();
        }
        this.socket.write(request);
    }

    private async getFreePort(): Promise<number> {
        return new Promise(resolve => {
            const ser = net.createServer();
            ser.listen(0, "127.0.0.1", () => {
                const addr = ser.address();
                ser.close();
                resolve((addr as net.AddressInfo).port);
            });
        });
    }

    private async delay(ms: number): Promise<unknown> {
        return new Promise((resolve) => setTimeout(resolve, ms));
    }

    private  spawnServer(args: string[]) {
        if (this.exeName === "Compose.exe") {
            process.env['PATH'] = path.join(this.omlRoot, '/hwx/bin/', this.platForm) + ";" +
                path.join(this.omlRoot, '/hw/bin/', this.platForm) + ";" + process.env['HW_FRAMEWORK'] +
                path.join('/bin/', this.platForm) + ";" + process.env['PATH'];
        }
        else {
            process.env['PATH'] = process.env['PATH'] + ";" +
                path.join(this.omlThirdParty, String(process.env['OML_INTEL_COMP'])) + ";" +
                path.join(this.omlThirdParty, String(process.env['OML_INTEL_MKL'])) + ";" +
                path.join(this.omlThirdParty, String(process.env['OML_FFTW'])) + ";" +
                path.join(this.omlThirdParty, String(process.env['OML_MATIO'])) + ";" +
                path.join(this.omlThirdParty, String(process.env['OML_HDF'])) + ";" +
                path.join(this.omlThirdParty, String(process.env['OML_QHULL']));
            console.log(process.env['PATH']);
        }
        let term =  this.createTerminal(args);
        return term?.processId;
    }

    private createTerminal(args: string[]): any {

        if (vscode.window.terminals.length > 0) {
            if (vscode.window.activeTerminal?.name === 'oml') {
                return vscode.window.activeTerminal;
            }
        }
        const ws = workspace.workspaceFolders;
        const wsPath = ws ? ws[0].uri.fsPath : undefined;
        const terminalOptions: TerminalOptions = {
            name: 'oml',
            shellPath: this.omlPath,
            shellArgs: args,
            cwd: wsPath,
            env: {
                ['PATH']: process.env['PATH'],
                ['QT_QPA_PLATFORM_PLUGIN_PATH']: process.env['QT_QPA_PLATFORM_PLUGIN_PATH'],
                ['OML_APPDIR']: process.env['OML_APPDIR'],
                ['OML_HELP']: process.env['OML_HELP'],
                ['ALTAIR_HOME']: process.env['ALTAIR_HOME'],
                ['HW_ROOTDIR']: process.env['HW_ROOTDIR'],
                ['HW_UNITY_ROOTDIR']: process.env['HW_UNITY_ROOTDIR'],
                ['HW_FRAMEWORK']: process.env['HW_FRAMEWORK']
            }
        };

        try {
            omlTerminal = window.createTerminal(terminalOptions);
            omlTerminal.show(false);
        } catch (error) {
            window.showErrorMessage("Could not create OML terminal.");
        }
        return omlTerminal;
    }

    public deleteTerminal(term: vscode.Terminal): void {
        if (term === omlTerminal) {
            omlTerminal.dispose();
            omlTerminal = undefined;
            this.isRunning = false;;
        }
    }

    public async runInTermianl(command: string) {
        if (omlTerminal) {
            omlTerminal.sendText(command);
            omlTerminal.show();
            await vscode.commands.executeCommand('workbench.action.terminal.scrollToBottom');
        }
    }

    private async startConnection(): Promise<LanguageClient> {
        async function didOpenTextDocument(document: TextDocument) {
            if (document.languageId !== 'oml' && document.languageId !== 'OML') {
                return;
            }
        }
        return await this.createClient();
    }

    public async restartConnection() {
        this.closeConnection();
        this.client = {} as LanguageClient;
        this.startConnection();
    }

    private closeConnection(): Thenable<void> {
        this.isRunning = false;
        const promises: Thenable<void>[] = [];
        promises.push(this.client.stop());
        return Promise.all(promises).then();
    }

    private async createClient(): Promise<LanguageClient> {

        this.outputChannel = window.createOutputChannel('OML');

        const port = await this.getFreePort();
        if (this.exeName === "Compose.exe") {
            const args = ['-withgui', '-toolbox', '-c', 'Compose', '-vscode-server', '-port', port.toString()];
             this.spawnServer(args);
        }
        else {
            const args = ['-vscode-server', '-port', port.toString()];
             this.spawnServer(args);
        }

        const serverOption = () => new Promise<StreamInfo>((resolve, reject) => {

            let socket = net.connect(port, 'localhost', () => {
                let result: StreamInfo = {
                    writer: socket,
                    reader: socket
                };
                console.log("Connection Established");
                return resolve(result);
            });
            this.socket = socket;
            socket.addListener("data", (chunk: Buffer) => this.recvData(chunk));
            socket.on('end', () => console.log('client disconnected from OML server'));
            socket.on('error', (e) => {
                console.log(e.message);
            });
        });

        const clientOptions: LanguageClientOptions = {
            // Register the server for OML documents
            documentSelector: [{ scheme: 'file', language: 'oml' }],
            outputChannel: this.outputChannel,
            synchronize: {
                // Notify the server about file changes to '.clientrc files contained in the workspace
                fileEvents: workspace.createFileSystemWatcher('**/.clientrc')
            }
        };
         // let the server warm up before client is created
        await this.delay(6000);
        let client = new LanguageClient('oml', 'OML Language Server', serverOption, clientOptions);
        extensionContext.subscriptions.push(client);
        this.setClient(client);
        try {
            await client.start();
        }
        catch (e) {
            vscode.window.showErrorMessage('Could not start the OML language server. Pls make sure the enviornment is set');
        }
        return client;
    }
}
