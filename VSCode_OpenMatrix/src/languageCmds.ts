import * as vscode from 'vscode';
import {langConn} from './extension';

export async function newFileDocument() {
    vscode.workspace.openTextDocument({language: 'oml'}).then((v) => vscode.window.showTextDocument(v));
}

export async function runActiveFile()
{
    if(! langConn.isConnectionAlive()) {
        langConn.restartConnection();
    }
    const actDoc = vscode.window.activeTextEditor?.document;
    if (!actDoc || actDoc.languageId !== 'oml') {
        return;
    }

    /*if (actDoc.isUntitled) {
        let text = actDoc.getText();
        let data: string = 'evalstring,' + text;
        langConn.sendRequest(data);
    }
    else {
        const isSaved = await actDoc.save();
        if (!isSaved) {
            vscode.window.showErrorMessage('Document could not be Saved. Pls check for write permission.');
        }
        let filePath: string = actDoc.uri.fsPath;
        let data: string = 'evalfile,' + filePath;
        langConn.sendRequest(data);
    }*/

    if (actDoc.isUntitled) 
    {
        let text = actDoc.getText();
        let data: string = 'evalstring,' + text;
        langConn.sendRequest(data);
        return;
    }
    await actDoc.save();
    let filePath: string = actDoc.uri.fsPath;
    let command: string = 'run(\'' + filePath + '\')';
    langConn.runInTermianl(command);
}

export async function stopExecution()
{
    if(langConn.isConnectionAlive()) {
        let data: string = 'interrupt,';
        langConn.sendRequest(data);
    }
}