OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6851299) q[0];
sx q[0];
rz(-0.87378341) q[0];
sx q[0];
rz(-1.1377347) q[0];
rz(-3.0955834) q[1];
sx q[1];
rz(0.53541056) q[1];
sx q[1];
rz(14.169197) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41661501) q[0];
sx q[0];
rz(-3.0862245) q[0];
sx q[0];
rz(2.2105818) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.12243226) q[2];
sx q[2];
rz(-0.85683595) q[2];
sx q[2];
rz(3.0453504) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5933758) q[1];
sx q[1];
rz(-2.0228099) q[1];
sx q[1];
rz(1.3424681) q[1];
x q[2];
rz(-3.0858062) q[3];
sx q[3];
rz(-1.7100289) q[3];
sx q[3];
rz(0.040797357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5753182) q[2];
sx q[2];
rz(-2.7988269) q[2];
sx q[2];
rz(2.9052367) q[2];
rz(-0.44860873) q[3];
sx q[3];
rz(-1.4834504) q[3];
sx q[3];
rz(2.1382704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3818632) q[0];
sx q[0];
rz(-2.6549082) q[0];
sx q[0];
rz(-0.78616649) q[0];
rz(-2.3568514) q[1];
sx q[1];
rz(-1.6678383) q[1];
sx q[1];
rz(3.0395708) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52962063) q[0];
sx q[0];
rz(-2.1054322) q[0];
sx q[0];
rz(-1.3360587) q[0];
rz(-0.78418248) q[2];
sx q[2];
rz(-1.3578556) q[2];
sx q[2];
rz(0.41894693) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1264915) q[1];
sx q[1];
rz(-0.41200629) q[1];
sx q[1];
rz(2.4088349) q[1];
rz(0.24436538) q[3];
sx q[3];
rz(-2.3373342) q[3];
sx q[3];
rz(1.3028631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.027779) q[2];
sx q[2];
rz(-1.6543829) q[2];
sx q[2];
rz(0.2629183) q[2];
rz(1.3024088) q[3];
sx q[3];
rz(-1.115256) q[3];
sx q[3];
rz(-2.3579679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1394434) q[0];
sx q[0];
rz(-0.0037010598) q[0];
sx q[0];
rz(-1.8078467) q[0];
rz(-1.3847146) q[1];
sx q[1];
rz(-0.92888558) q[1];
sx q[1];
rz(0.76470107) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0686808) q[0];
sx q[0];
rz(-0.31038767) q[0];
sx q[0];
rz(-2.458802) q[0];
rz(-2.7638859) q[2];
sx q[2];
rz(-1.3481082) q[2];
sx q[2];
rz(-1.6645704) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0460658) q[1];
sx q[1];
rz(-2.8316759) q[1];
sx q[1];
rz(-0.97280963) q[1];
x q[2];
rz(-2.2883203) q[3];
sx q[3];
rz(-1.3382592) q[3];
sx q[3];
rz(-2.5465089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9646405) q[2];
sx q[2];
rz(-1.8884337) q[2];
sx q[2];
rz(-2.9580252) q[2];
rz(-0.14487264) q[3];
sx q[3];
rz(-1.262007) q[3];
sx q[3];
rz(-0.22411331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1862828) q[0];
sx q[0];
rz(-1.1306385) q[0];
sx q[0];
rz(2.8593707) q[0];
rz(1.9918293) q[1];
sx q[1];
rz(-1.3590004) q[1];
sx q[1];
rz(1.6414292) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1896797) q[0];
sx q[0];
rz(-1.3985347) q[0];
sx q[0];
rz(-2.7776615) q[0];
x q[1];
rz(2.8871782) q[2];
sx q[2];
rz(-1.0378278) q[2];
sx q[2];
rz(2.0526469) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.73622656) q[1];
sx q[1];
rz(-0.42441503) q[1];
sx q[1];
rz(-2.2176803) q[1];
x q[2];
rz(0.30434609) q[3];
sx q[3];
rz(-1.4831717) q[3];
sx q[3];
rz(-2.9612142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.942261) q[2];
sx q[2];
rz(-1.4456238) q[2];
sx q[2];
rz(1.4873827) q[2];
rz(-2.2516294) q[3];
sx q[3];
rz(-1.7421937) q[3];
sx q[3];
rz(-0.14931211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99059659) q[0];
sx q[0];
rz(-2.6166333) q[0];
sx q[0];
rz(-1.6816444) q[0];
rz(1.2851985) q[1];
sx q[1];
rz(-1.3672914) q[1];
sx q[1];
rz(-2.6885414) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.995034) q[0];
sx q[0];
rz(-1.5164598) q[0];
sx q[0];
rz(-1.2666235) q[0];
x q[1];
rz(-1.4905246) q[2];
sx q[2];
rz(-0.58387127) q[2];
sx q[2];
rz(-0.50498) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.510524) q[1];
sx q[1];
rz(-1.9434526) q[1];
sx q[1];
rz(-1.371553) q[1];
rz(-pi) q[2];
rz(1.1207709) q[3];
sx q[3];
rz(-1.1965568) q[3];
sx q[3];
rz(0.14412021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4840661) q[2];
sx q[2];
rz(-2.6675197) q[2];
sx q[2];
rz(-1.5849812) q[2];
rz(1.5629684) q[3];
sx q[3];
rz(-1.2789187) q[3];
sx q[3];
rz(2.1141619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9180561) q[0];
sx q[0];
rz(-2.1503088) q[0];
sx q[0];
rz(1.1728485) q[0];
rz(2.6808443) q[1];
sx q[1];
rz(-1.2604424) q[1];
sx q[1];
rz(-1.4985098) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39678412) q[0];
sx q[0];
rz(-0.61273324) q[0];
sx q[0];
rz(-1.5694797) q[0];
x q[1];
rz(0.66835673) q[2];
sx q[2];
rz(-0.41497013) q[2];
sx q[2];
rz(0.1842095) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9869796) q[1];
sx q[1];
rz(-1.8334532) q[1];
sx q[1];
rz(-3.0820552) q[1];
rz(1.40703) q[3];
sx q[3];
rz(-2.4330045) q[3];
sx q[3];
rz(2.6276789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1341165) q[2];
sx q[2];
rz(-1.7191929) q[2];
sx q[2];
rz(-2.9134992) q[2];
rz(-0.54276931) q[3];
sx q[3];
rz(-2.1398862) q[3];
sx q[3];
rz(-0.91226474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0246564) q[0];
sx q[0];
rz(-1.0448562) q[0];
sx q[0];
rz(-0.60633099) q[0];
rz(1.2888651) q[1];
sx q[1];
rz(-1.5325129) q[1];
sx q[1];
rz(0.85561633) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.037362075) q[0];
sx q[0];
rz(-2.3920569) q[0];
sx q[0];
rz(-0.61893344) q[0];
rz(0.77328859) q[2];
sx q[2];
rz(-1.5097364) q[2];
sx q[2];
rz(1.9715015) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.85904494) q[1];
sx q[1];
rz(-0.63524109) q[1];
sx q[1];
rz(2.7905491) q[1];
x q[2];
rz(-1.0904013) q[3];
sx q[3];
rz(-1.3452936) q[3];
sx q[3];
rz(0.71798872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3768846) q[2];
sx q[2];
rz(-2.7089684) q[2];
sx q[2];
rz(0.077733668) q[2];
rz(-3.0149095) q[3];
sx q[3];
rz(-2.099791) q[3];
sx q[3];
rz(-2.3104987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7938101) q[0];
sx q[0];
rz(-0.3599444) q[0];
sx q[0];
rz(0.98156324) q[0];
rz(1.7836102) q[1];
sx q[1];
rz(-0.49109083) q[1];
sx q[1];
rz(0.29409274) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9413853) q[0];
sx q[0];
rz(-1.6305171) q[0];
sx q[0];
rz(-2.8325547) q[0];
rz(-pi) q[1];
rz(0.63801672) q[2];
sx q[2];
rz(-2.6891481) q[2];
sx q[2];
rz(2.5663301) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24178019) q[1];
sx q[1];
rz(-2.2072148) q[1];
sx q[1];
rz(2.7905725) q[1];
rz(-1.4979248) q[3];
sx q[3];
rz(-1.6207063) q[3];
sx q[3];
rz(-2.2467256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.10245094) q[2];
sx q[2];
rz(-2.4852018) q[2];
sx q[2];
rz(1.0605109) q[2];
rz(0.55824009) q[3];
sx q[3];
rz(-2.5679936) q[3];
sx q[3];
rz(3.1225045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3625665) q[0];
sx q[0];
rz(-2.9988852) q[0];
sx q[0];
rz(2.3574164) q[0];
rz(1.7991637) q[1];
sx q[1];
rz(-1.289184) q[1];
sx q[1];
rz(0.69650355) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9331037) q[0];
sx q[0];
rz(-1.7835) q[0];
sx q[0];
rz(-1.8038294) q[0];
rz(2.0663459) q[2];
sx q[2];
rz(-2.181567) q[2];
sx q[2];
rz(-1.3411759) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.15882684) q[1];
sx q[1];
rz(-2.4707795) q[1];
sx q[1];
rz(2.4300086) q[1];
x q[2];
rz(0.83902766) q[3];
sx q[3];
rz(-0.19908842) q[3];
sx q[3];
rz(1.3092878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8037618) q[2];
sx q[2];
rz(-0.75118128) q[2];
sx q[2];
rz(1.8776228) q[2];
rz(1.3761282) q[3];
sx q[3];
rz(-2.2270146) q[3];
sx q[3];
rz(-1.5862563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2176168) q[0];
sx q[0];
rz(-0.042984977) q[0];
sx q[0];
rz(1.6643583) q[0];
rz(-2.2337275) q[1];
sx q[1];
rz(-1.5217179) q[1];
sx q[1];
rz(0.97000617) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5610301) q[0];
sx q[0];
rz(-0.99107689) q[0];
sx q[0];
rz(0.28663488) q[0];
rz(0.97899784) q[2];
sx q[2];
rz(-1.819397) q[2];
sx q[2];
rz(-0.73038855) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2737588) q[1];
sx q[1];
rz(-0.9504488) q[1];
sx q[1];
rz(1.7628487) q[1];
rz(-pi) q[2];
rz(-0.72204263) q[3];
sx q[3];
rz(-1.0580499) q[3];
sx q[3];
rz(-1.7554612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.19405356) q[2];
sx q[2];
rz(-2.2668362) q[2];
sx q[2];
rz(-0.27210316) q[2];
rz(2.2943606) q[3];
sx q[3];
rz(-0.50744414) q[3];
sx q[3];
rz(1.0890755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.465441) q[0];
sx q[0];
rz(-2.0068824) q[0];
sx q[0];
rz(-1.6750499) q[0];
rz(-0.33879406) q[1];
sx q[1];
rz(-1.6962961) q[1];
sx q[1];
rz(1.2512339) q[1];
rz(1.2180381) q[2];
sx q[2];
rz(-2.259476) q[2];
sx q[2];
rz(-1.8720575) q[2];
rz(1.1477039) q[3];
sx q[3];
rz(-1.4948899) q[3];
sx q[3];
rz(-0.75567452) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
