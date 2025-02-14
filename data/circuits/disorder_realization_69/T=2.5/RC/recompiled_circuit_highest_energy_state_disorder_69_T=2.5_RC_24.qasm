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
rz(1.5975098) q[0];
sx q[0];
rz(-1.37356) q[0];
sx q[0];
rz(-2.1231667) q[0];
rz(2.6234558) q[1];
sx q[1];
rz(-0.75704804) q[1];
sx q[1];
rz(2.5098324) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9874532) q[0];
sx q[0];
rz(-1.69603) q[0];
sx q[0];
rz(1.7787329) q[0];
rz(1.2680797) q[2];
sx q[2];
rz(-1.5945425) q[2];
sx q[2];
rz(1.0645006) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.85447299) q[1];
sx q[1];
rz(-1.8504256) q[1];
sx q[1];
rz(-0.92052144) q[1];
rz(0.092436304) q[3];
sx q[3];
rz(-1.2458003) q[3];
sx q[3];
rz(0.50918885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.53039256) q[2];
sx q[2];
rz(-2.7860614) q[2];
sx q[2];
rz(-1.6543039) q[2];
rz(-2.2330331) q[3];
sx q[3];
rz(-0.23659758) q[3];
sx q[3];
rz(-1.0049741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0205883) q[0];
sx q[0];
rz(-1.7089184) q[0];
sx q[0];
rz(0.16227907) q[0];
rz(0.24457112) q[1];
sx q[1];
rz(-1.9385612) q[1];
sx q[1];
rz(-2.8499106) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2431716) q[0];
sx q[0];
rz(-2.7968639) q[0];
sx q[0];
rz(-2.5401462) q[0];
rz(-2.2333916) q[2];
sx q[2];
rz(-2.78763) q[2];
sx q[2];
rz(2.5733054) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.32714601) q[1];
sx q[1];
rz(-1.3400048) q[1];
sx q[1];
rz(-0.77951851) q[1];
x q[2];
rz(-1.9883745) q[3];
sx q[3];
rz(-1.1912701) q[3];
sx q[3];
rz(-2.1714927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.67819277) q[2];
sx q[2];
rz(-2.2726111) q[2];
sx q[2];
rz(2.0657067) q[2];
rz(1.2470657) q[3];
sx q[3];
rz(-1.5445292) q[3];
sx q[3];
rz(1.4472848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4397044) q[0];
sx q[0];
rz(-2.0151558) q[0];
sx q[0];
rz(0.18371789) q[0];
rz(1.6366929) q[1];
sx q[1];
rz(-0.74438649) q[1];
sx q[1];
rz(2.9768129) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67023428) q[0];
sx q[0];
rz(-0.58588282) q[0];
sx q[0];
rz(-1.5536932) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7585331) q[2];
sx q[2];
rz(-0.7128517) q[2];
sx q[2];
rz(1.0733611) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7565733) q[1];
sx q[1];
rz(-2.1655271) q[1];
sx q[1];
rz(-1.1043332) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.80071) q[3];
sx q[3];
rz(-1.0629144) q[3];
sx q[3];
rz(-1.1424049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.1068716) q[2];
sx q[2];
rz(-1.9811337) q[2];
sx q[2];
rz(-1.1064233) q[2];
rz(2.4404081) q[3];
sx q[3];
rz(-1.3283575) q[3];
sx q[3];
rz(-0.4655233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.768854) q[0];
sx q[0];
rz(-1.1356069) q[0];
sx q[0];
rz(-2.1970774) q[0];
rz(1.6230029) q[1];
sx q[1];
rz(-0.98097643) q[1];
sx q[1];
rz(-1.5400344) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0991411) q[0];
sx q[0];
rz(-0.97223982) q[0];
sx q[0];
rz(0.70229806) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3643199) q[2];
sx q[2];
rz(-0.16675719) q[2];
sx q[2];
rz(0.47698944) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9303927) q[1];
sx q[1];
rz(-2.2477149) q[1];
sx q[1];
rz(-2.6571214) q[1];
rz(-1.6155195) q[3];
sx q[3];
rz(-2.3207242) q[3];
sx q[3];
rz(-2.5673411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.010178415) q[2];
sx q[2];
rz(-2.5203036) q[2];
sx q[2];
rz(2.9739001) q[2];
rz(0.0532648) q[3];
sx q[3];
rz(-2.0361418) q[3];
sx q[3];
rz(1.7648599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5662956) q[0];
sx q[0];
rz(-3.0360041) q[0];
sx q[0];
rz(2.8422624) q[0];
rz(-2.7686367) q[1];
sx q[1];
rz(-1.2495709) q[1];
sx q[1];
rz(-1.6711055) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9540967) q[0];
sx q[0];
rz(-0.72802714) q[0];
sx q[0];
rz(-2.4690829) q[0];
rz(-pi) q[1];
rz(0.23678933) q[2];
sx q[2];
rz(-2.6663187) q[2];
sx q[2];
rz(-1.4314112) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.336323) q[1];
sx q[1];
rz(-2.0741558) q[1];
sx q[1];
rz(0.663228) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8138795) q[3];
sx q[3];
rz(-0.63589261) q[3];
sx q[3];
rz(2.2896374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2028929) q[2];
sx q[2];
rz(-2.4788269) q[2];
sx q[2];
rz(2.0311484) q[2];
rz(-0.86197305) q[3];
sx q[3];
rz(-1.6317261) q[3];
sx q[3];
rz(-1.2601669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92457572) q[0];
sx q[0];
rz(-0.68704263) q[0];
sx q[0];
rz(-0.4050912) q[0];
rz(-2.8054667) q[1];
sx q[1];
rz(-1.4210217) q[1];
sx q[1];
rz(0.95132557) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66621214) q[0];
sx q[0];
rz(-2.6586091) q[0];
sx q[0];
rz(0.57421143) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4053965) q[2];
sx q[2];
rz(-1.412782) q[2];
sx q[2];
rz(-2.406213) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3094294) q[1];
sx q[1];
rz(-1.0069261) q[1];
sx q[1];
rz(1.840074) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5441465) q[3];
sx q[3];
rz(-0.33388147) q[3];
sx q[3];
rz(0.88874879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.090791) q[2];
sx q[2];
rz(-1.4046706) q[2];
sx q[2];
rz(0.23078272) q[2];
rz(2.9595621) q[3];
sx q[3];
rz(-2.5735276) q[3];
sx q[3];
rz(-0.52869421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35293216) q[0];
sx q[0];
rz(-3.1020628) q[0];
sx q[0];
rz(0.93210644) q[0];
rz(-0.63356361) q[1];
sx q[1];
rz(-2.0220058) q[1];
sx q[1];
rz(1.8036141) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.103823) q[0];
sx q[0];
rz(-1.4892206) q[0];
sx q[0];
rz(-1.6768811) q[0];
x q[1];
rz(-1.6858844) q[2];
sx q[2];
rz(-1.8427927) q[2];
sx q[2];
rz(0.65036648) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4949957) q[1];
sx q[1];
rz(-0.75237583) q[1];
sx q[1];
rz(-0.97480358) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0831117) q[3];
sx q[3];
rz(-0.49904682) q[3];
sx q[3];
rz(2.8935695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.46447095) q[2];
sx q[2];
rz(-1.2268343) q[2];
sx q[2];
rz(-1.9390437) q[2];
rz(2.7770212) q[3];
sx q[3];
rz(-2.4489844) q[3];
sx q[3];
rz(0.83474368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6637591) q[0];
sx q[0];
rz(-0.57357016) q[0];
sx q[0];
rz(0.17663503) q[0];
rz(-2.7015576) q[1];
sx q[1];
rz(-1.9562079) q[1];
sx q[1];
rz(-2.1766591) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4266708) q[0];
sx q[0];
rz(-1.3765125) q[0];
sx q[0];
rz(1.4358836) q[0];
rz(-pi) q[1];
rz(2.1148483) q[2];
sx q[2];
rz(-2.1026582) q[2];
sx q[2];
rz(3.056125) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7495887) q[1];
sx q[1];
rz(-0.67396213) q[1];
sx q[1];
rz(0.96629179) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3163525) q[3];
sx q[3];
rz(-1.644205) q[3];
sx q[3];
rz(1.0496828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9289916) q[2];
sx q[2];
rz(-1.6182199) q[2];
sx q[2];
rz(3.1316481) q[2];
rz(-3.0250004) q[3];
sx q[3];
rz(-0.3392342) q[3];
sx q[3];
rz(0.54178437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56124878) q[0];
sx q[0];
rz(-1.9191701) q[0];
sx q[0];
rz(0.62193459) q[0];
rz(1.4382582) q[1];
sx q[1];
rz(-2.56918) q[1];
sx q[1];
rz(-2.4780746) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2885482) q[0];
sx q[0];
rz(-0.85381928) q[0];
sx q[0];
rz(-1.8621481) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2873257) q[2];
sx q[2];
rz(-1.4525692) q[2];
sx q[2];
rz(-2.9663939) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.89412824) q[1];
sx q[1];
rz(-0.93527764) q[1];
sx q[1];
rz(0.39860098) q[1];
rz(-pi) q[2];
rz(-2.8264753) q[3];
sx q[3];
rz(-2.8459918) q[3];
sx q[3];
rz(0.35741266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5180987) q[2];
sx q[2];
rz(-1.3891209) q[2];
sx q[2];
rz(0.36250472) q[2];
rz(2.6643961) q[3];
sx q[3];
rz(-2.0878744) q[3];
sx q[3];
rz(1.1227192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0839194) q[0];
sx q[0];
rz(-0.73129439) q[0];
sx q[0];
rz(-0.90173632) q[0];
rz(-2.4665191) q[1];
sx q[1];
rz(-0.98552862) q[1];
sx q[1];
rz(0.14428446) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1146204) q[0];
sx q[0];
rz(-1.2204224) q[0];
sx q[0];
rz(3.0536985) q[0];
x q[1];
rz(0.4298019) q[2];
sx q[2];
rz(-0.24082213) q[2];
sx q[2];
rz(-2.4483829) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.212483) q[1];
sx q[1];
rz(-1.3204201) q[1];
sx q[1];
rz(-0.65156619) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7561109) q[3];
sx q[3];
rz(-0.47244888) q[3];
sx q[3];
rz(-0.89490283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5054063) q[2];
sx q[2];
rz(-1.212684) q[2];
sx q[2];
rz(0.20067659) q[2];
rz(1.1096654) q[3];
sx q[3];
rz(-2.6789013) q[3];
sx q[3];
rz(-0.5717352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5526445) q[0];
sx q[0];
rz(-2.1727967) q[0];
sx q[0];
rz(-1.1217242) q[0];
rz(-2.6523392) q[1];
sx q[1];
rz(-1.4514634) q[1];
sx q[1];
rz(-1.0101752) q[1];
rz(1.1560925) q[2];
sx q[2];
rz(-1.4857875) q[2];
sx q[2];
rz(-0.93405741) q[2];
rz(-0.060754178) q[3];
sx q[3];
rz(-2.698632) q[3];
sx q[3];
rz(-2.5198577) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
