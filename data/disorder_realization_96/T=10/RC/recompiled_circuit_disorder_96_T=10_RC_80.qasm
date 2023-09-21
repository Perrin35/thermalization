OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.049790073) q[0];
sx q[0];
rz(-0.12806211) q[0];
sx q[0];
rz(0.81737104) q[0];
rz(0.983239) q[1];
sx q[1];
rz(2.6020738) q[1];
sx q[1];
rz(10.625216) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96785986) q[0];
sx q[0];
rz(-1.1636486) q[0];
sx q[0];
rz(-1.3074387) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3517411) q[2];
sx q[2];
rz(-0.54974216) q[2];
sx q[2];
rz(1.4866536) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.22908902) q[1];
sx q[1];
rz(-0.87516811) q[1];
sx q[1];
rz(1.9181262) q[1];
rz(0.79521631) q[3];
sx q[3];
rz(-1.3391558) q[3];
sx q[3];
rz(-2.8455646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1203221) q[2];
sx q[2];
rz(-0.56818429) q[2];
sx q[2];
rz(-1.5585287) q[2];
rz(0.99672404) q[3];
sx q[3];
rz(-0.45209) q[3];
sx q[3];
rz(-2.7157917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5581756) q[0];
sx q[0];
rz(-0.75220627) q[0];
sx q[0];
rz(-3.0875207) q[0];
rz(1.1955098) q[1];
sx q[1];
rz(-2.1046488) q[1];
sx q[1];
rz(-0.53584677) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5553404) q[0];
sx q[0];
rz(-0.99389168) q[0];
sx q[0];
rz(0.14848407) q[0];
rz(-pi) q[1];
rz(-0.29157721) q[2];
sx q[2];
rz(-2.4107286) q[2];
sx q[2];
rz(2.4183194) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9745001) q[1];
sx q[1];
rz(-1.316861) q[1];
sx q[1];
rz(-0.42029917) q[1];
rz(3.1090118) q[3];
sx q[3];
rz(-1.7136646) q[3];
sx q[3];
rz(-2.8829765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0216996) q[2];
sx q[2];
rz(-1.0571486) q[2];
sx q[2];
rz(-2.1832441) q[2];
rz(3.0751394) q[3];
sx q[3];
rz(-1.5575912) q[3];
sx q[3];
rz(0.44979969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41985837) q[0];
sx q[0];
rz(-1.9323213) q[0];
sx q[0];
rz(-0.15047519) q[0];
rz(-0.45723215) q[1];
sx q[1];
rz(-0.91232863) q[1];
sx q[1];
rz(3.1157852) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8919825) q[0];
sx q[0];
rz(-1.550807) q[0];
sx q[0];
rz(1.32094) q[0];
rz(-pi) q[1];
rz(3.0518968) q[2];
sx q[2];
rz(-1.8811474) q[2];
sx q[2];
rz(2.2657564) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3734696) q[1];
sx q[1];
rz(-0.79263955) q[1];
sx q[1];
rz(-2.5440689) q[1];
x q[2];
rz(-1.8886391) q[3];
sx q[3];
rz(-2.7405973) q[3];
sx q[3];
rz(-2.3272115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1228483) q[2];
sx q[2];
rz(-2.7763425) q[2];
sx q[2];
rz(-2.5562111) q[2];
rz(0.18150005) q[3];
sx q[3];
rz(-1.8094962) q[3];
sx q[3];
rz(-1.5766778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.240775) q[0];
sx q[0];
rz(-2.5188991) q[0];
sx q[0];
rz(2.9649819) q[0];
rz(-0.88090849) q[1];
sx q[1];
rz(-1.0662339) q[1];
sx q[1];
rz(2.6054629) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6648383) q[0];
sx q[0];
rz(-1.9739445) q[0];
sx q[0];
rz(2.320015) q[0];
x q[1];
rz(-1.6279814) q[2];
sx q[2];
rz(-1.3344889) q[2];
sx q[2];
rz(-1.1255217) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9163461) q[1];
sx q[1];
rz(-0.9496453) q[1];
sx q[1];
rz(-2.0398554) q[1];
x q[2];
rz(2.4651277) q[3];
sx q[3];
rz(-1.3389265) q[3];
sx q[3];
rz(-1.3456618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6716016) q[2];
sx q[2];
rz(-1.4208379) q[2];
sx q[2];
rz(-1.0446576) q[2];
rz(0.70703834) q[3];
sx q[3];
rz(-2.1576594) q[3];
sx q[3];
rz(0.35693359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2351284) q[0];
sx q[0];
rz(-1.8179853) q[0];
sx q[0];
rz(2.0902324) q[0];
rz(-1.6479187) q[1];
sx q[1];
rz(-2.5525679) q[1];
sx q[1];
rz(0.043118127) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9263822) q[0];
sx q[0];
rz(-0.97391093) q[0];
sx q[0];
rz(-2.9125288) q[0];
x q[1];
rz(-2.0448858) q[2];
sx q[2];
rz(-2.5708963) q[2];
sx q[2];
rz(1.8622461) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9777898) q[1];
sx q[1];
rz(-1.3063352) q[1];
sx q[1];
rz(-2.8603641) q[1];
rz(-pi) q[2];
rz(2.352036) q[3];
sx q[3];
rz(-1.6766747) q[3];
sx q[3];
rz(-1.7136128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0126426) q[2];
sx q[2];
rz(-2.6439715) q[2];
sx q[2];
rz(-2.7894003) q[2];
rz(-2.5514065) q[3];
sx q[3];
rz(-2.6679109) q[3];
sx q[3];
rz(-0.56110704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6234289) q[0];
sx q[0];
rz(-1.9264899) q[0];
sx q[0];
rz(1.9859001) q[0];
rz(-2.39134) q[1];
sx q[1];
rz(-2.2022088) q[1];
sx q[1];
rz(-1.0587143) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6569865) q[0];
sx q[0];
rz(-1.3625506) q[0];
sx q[0];
rz(2.7013742) q[0];
x q[1];
rz(-2.6201453) q[2];
sx q[2];
rz(-1.2505184) q[2];
sx q[2];
rz(-2.5495868) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2211654) q[1];
sx q[1];
rz(-1.3472918) q[1];
sx q[1];
rz(1.8799053) q[1];
x q[2];
rz(-0.70763564) q[3];
sx q[3];
rz(-1.71873) q[3];
sx q[3];
rz(1.9572452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6713312) q[2];
sx q[2];
rz(-1.7241717) q[2];
sx q[2];
rz(-1.2188101) q[2];
rz(1.1550711) q[3];
sx q[3];
rz(-0.23854908) q[3];
sx q[3];
rz(-1.7003805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041615151) q[0];
sx q[0];
rz(-1.8247373) q[0];
sx q[0];
rz(1.51145) q[0];
rz(1.7639683) q[1];
sx q[1];
rz(-2.8306077) q[1];
sx q[1];
rz(-2.2999433) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8182897) q[0];
sx q[0];
rz(-1.0536195) q[0];
sx q[0];
rz(1.3099758) q[0];
rz(1.7633348) q[2];
sx q[2];
rz(-2.1888869) q[2];
sx q[2];
rz(-0.75380177) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.48982606) q[1];
sx q[1];
rz(-0.090432743) q[1];
sx q[1];
rz(-1.0033146) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3326725) q[3];
sx q[3];
rz(-2.3861285) q[3];
sx q[3];
rz(-1.8734224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1402309) q[2];
sx q[2];
rz(-1.3683687) q[2];
sx q[2];
rz(-3.0916396) q[2];
rz(2.4800381) q[3];
sx q[3];
rz(-0.52246061) q[3];
sx q[3];
rz(0.18930999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2426303) q[0];
sx q[0];
rz(-2.1147418) q[0];
sx q[0];
rz(1.4021953) q[0];
rz(3.0461123) q[1];
sx q[1];
rz(-1.9752558) q[1];
sx q[1];
rz(-0.41762525) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1091052) q[0];
sx q[0];
rz(-2.0628477) q[0];
sx q[0];
rz(-2.0672654) q[0];
x q[1];
rz(1.8629486) q[2];
sx q[2];
rz(-2.7610965) q[2];
sx q[2];
rz(-0.7064864) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0954674) q[1];
sx q[1];
rz(-2.1324665) q[1];
sx q[1];
rz(-2.5820877) q[1];
rz(-pi) q[2];
rz(-0.29406677) q[3];
sx q[3];
rz(-1.9293647) q[3];
sx q[3];
rz(-0.29953526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3461356) q[2];
sx q[2];
rz(-1.6980349) q[2];
sx q[2];
rz(-1.0428838) q[2];
rz(0.67388326) q[3];
sx q[3];
rz(-1.4833114) q[3];
sx q[3];
rz(-0.90014443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3089356) q[0];
sx q[0];
rz(-1.7333663) q[0];
sx q[0];
rz(2.5119264) q[0];
rz(0.57485238) q[1];
sx q[1];
rz(-1.8300627) q[1];
sx q[1];
rz(2.1946857) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1614721) q[0];
sx q[0];
rz(-1.8959909) q[0];
sx q[0];
rz(-1.8878216) q[0];
x q[1];
rz(-0.097054585) q[2];
sx q[2];
rz(-1.5311054) q[2];
sx q[2];
rz(1.5904215) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7864101) q[1];
sx q[1];
rz(-1.2372412) q[1];
sx q[1];
rz(1.6598318) q[1];
rz(0.19091786) q[3];
sx q[3];
rz(-2.0058161) q[3];
sx q[3];
rz(0.85489475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5809014) q[2];
sx q[2];
rz(-1.1297444) q[2];
sx q[2];
rz(-1.8927195) q[2];
rz(-2.4272264) q[3];
sx q[3];
rz(-1.2759821) q[3];
sx q[3];
rz(-2.8760288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0666075) q[0];
sx q[0];
rz(-2.5191436) q[0];
sx q[0];
rz(0.65504909) q[0];
rz(-2.24522) q[1];
sx q[1];
rz(-2.2348149) q[1];
sx q[1];
rz(0.64430976) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33675942) q[0];
sx q[0];
rz(-1.5157489) q[0];
sx q[0];
rz(0.17382646) q[0];
rz(-pi) q[1];
rz(-2.2868025) q[2];
sx q[2];
rz(-1.3959179) q[2];
sx q[2];
rz(1.0027494) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2584553) q[1];
sx q[1];
rz(-1.2346134) q[1];
sx q[1];
rz(-1.8717481) q[1];
x q[2];
rz(0.44390042) q[3];
sx q[3];
rz(-1.1392986) q[3];
sx q[3];
rz(-0.75507009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7908988) q[2];
sx q[2];
rz(-1.472298) q[2];
sx q[2];
rz(-2.541686) q[2];
rz(-0.89896262) q[3];
sx q[3];
rz(-0.18342429) q[3];
sx q[3];
rz(-1.775734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8469289) q[0];
sx q[0];
rz(-1.3273205) q[0];
sx q[0];
rz(2.474665) q[0];
rz(-2.9121493) q[1];
sx q[1];
rz(-0.89090092) q[1];
sx q[1];
rz(0.13577239) q[1];
rz(2.1951998) q[2];
sx q[2];
rz(-1.6740435) q[2];
sx q[2];
rz(-1.2798535) q[2];
rz(-0.89623981) q[3];
sx q[3];
rz(-1.8825718) q[3];
sx q[3];
rz(2.752302) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];