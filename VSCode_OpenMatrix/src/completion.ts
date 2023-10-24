import * as vscode from 'vscode';
import {langConn} from './extension';

export class CompletionData {
    private keywords_:string[];
    private functions_:string[];

    constructor(first:string, last:string)
    {
        this.keywords_=[];
        this.functions_=[];
    }

    setKeywords(keywords:string) {
        this.keywords_ = keywords.split(' ');
    }
    setFunctions(functions:string) {
        this.functions_ = functions.split(' ');
    }
    getCompletionItems(target:string):string[]  {
        let regExp = RegExp('^' + target + '(?=).*');
        var completionList = this.keywords_.concat(this.functions_);
        let items = completionList.filter((str) => {
            return regExp.test(str);
        });
        items.sort((lhs, rhs) => {
            return lhs < rhs ? -1 : lhs > rhs ? 1 : 0;
        });
        return items;
    }
}
export class CompletionItemProvider implements vscode.CompletionItemProvider {

    private completionData:CompletionData;

    constructor(compData:CompletionData)
    {
        this.completionData = compData;
    }

    provideCompletionItems(
        document: vscode.TextDocument,
        position: vscode.Position,
        token: vscode.CancellationToken,
        completionContext: vscode.CompletionContext
    ): vscode.CompletionItem[] {
        
        const items: vscode.CompletionItem[] = [];
        if (token.isCancellationRequested) {
            return items;
        }
        const targetPosition = new vscode.Position(position.line, position.character - 1);
        const targetRange = document.getWordRangeAtPosition(targetPosition);
        const target = document.getText(targetRange);
        if(target.length > 1) {
            const matches = this.completionData.getCompletionItems(target);
            matches.forEach( (match) => {
                items.push(new vscode.CompletionItem(match, vscode.CompletionItemKind.Function));
            });
        }
        return items;
    }
}

export class SignatuareHelpProvider implements vscode.SignatureHelpProvider {

    private target_:string = '';
    private commas_:number = 0;

    async provideSignatureHelp(
        document: vscode.TextDocument,
        position: vscode.Position,
        token: vscode.CancellationToken,
        context: vscode.SignatureHelpContext
    ): Promise<vscode.SignatureHelp> {

        this.commas_ = 0;
        this.target_ = '';

        if (token.isCancellationRequested || context === undefined) {
            return new vscode.SignatureHelp();
        }

        if(context.triggerCharacter === '(') {
            const targetPosition = new vscode.Position(position.line, position.character - 1);
            const targetRange = document.getWordRangeAtPosition(targetPosition);
            this.target_ = document.getText(targetRange);
        }
        else  {
            // scan the line for target
            let lineText = document.lineAt(position.line).text;
            for(let i = position.character - 1; i > 0; --i ) {
                if(lineText[i] === '(') {
                    const targetPos = new vscode.Position(position.line, i-1);
                    const targetRange = document.getWordRangeAtPosition(targetPos);
                    this.target_ = document.getText(targetRange);
                    break;
                }
                else if(lineText[i] === ')') {
                    this.target_ = '';
                    break;
                }
                else if(lineText[i] === ',') {
                    ++this.commas_;
                }
            }
        }
        
        if(this.target_.length === 0) {
            return new vscode.SignatureHelp();
        }

        let data: string = 'signature,' + this.target_;
        await langConn.sendRequest(data);

        return new Promise<vscode.SignatureHelp>(resolve => {
            function provideSign(sign:string, commas:number) {
                let sigHelp = new vscode.SignatureHelp;
                sigHelp.activeSignature = 0;
                sigHelp.activeParameter = 0;
                let signs = sign.split(/(?<=\.)\n/);
                signs.forEach((data) => { 
                    let sinfo = new vscode.SignatureInformation(data);
                    let params;
                    if((params = data.match(/\((.*)\)/)) !== null) {
                        for(let i = 1; i < params.length; ++i) {
                            let param = params[i].split(',');
                            param.forEach(val => {
                                val.trim();
                                if(val.length) {
                                    sinfo.parameters.push(new vscode.ParameterInformation(val));
                                }
                            });
                            if(commas > 0 &&  commas === sinfo.parameters.length-1) {
                                sigHelp.activeSignature = sigHelp.signatures.length;
                                sigHelp.activeParameter = commas;
                            }
                        }
                    }
                    sigHelp.signatures.push(sinfo);
                });
                resolve(sigHelp);
            }

            langConn.addListener('signature', (sign:string) => {
                provideSign(sign, this.commas_);
            });
        });
    }
}

//export class OmlLanguageConfiguration implements vscode.LanguageConfiguration {}

/*export class SemanticTokensProvider implements vscode.DocumentSemanticTokensProvider {
    async provideDocumentSemanticTokens(document: vscode.TextDocument, token: vscode.CancellationToken): Promise<vscode.SemanticTokens> {
        const builder = new vscode.SemanticTokensBuilder();

    }
}*/