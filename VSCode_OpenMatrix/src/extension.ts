// The module 'vscode' contains the VS Code extensibility API
// Import the module and reference it with the alias vscode in your code below
import * as vscode from 'vscode';
import * as languageCmds from './languageCmds';
import { LanguageConnection } from "./languageConnection";

export let extensionContext:vscode.ExtensionContext;
export let langConn: LanguageConnection;

// This method is called when OML extension is activated
// The extension is activated the very first time the command is executed
export async function activate(context: vscode.ExtensionContext) {

	 // register commands from package.json
	 context.subscriptions.push(vscode.commands.registerCommand('oml.Run', languageCmds.runActiveFile));
	 context.subscriptions.push(vscode.commands.registerCommand('oml.Stop', languageCmds.stopExecution));
	 context.subscriptions.push(vscode.commands.registerCommand('oml.newFileDocument',languageCmds.newFileDocument));


	// Create the language connection and start the client.
	extensionContext = context;
	langConn = new LanguageConnection();
	context.subscriptions.push(langConn);
}

// This method is called when your extension is deactivated
export async function deactivate() 
{
	await langConn.dispose();
}
