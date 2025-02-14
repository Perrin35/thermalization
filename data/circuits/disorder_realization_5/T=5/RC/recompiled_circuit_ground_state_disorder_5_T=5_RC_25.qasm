OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.65052819) q[0];
sx q[0];
rz(-1.0125546) q[0];
sx q[0];
rz(-2.2192686) q[0];
rz(2.2946279) q[1];
sx q[1];
rz(-1.4743409) q[1];
sx q[1];
rz(2.9234731) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6309752) q[0];
sx q[0];
rz(-1.4790863) q[0];
sx q[0];
rz(2.6789078) q[0];
x q[1];
rz(-3.0384205) q[2];
sx q[2];
rz(-2.7617447) q[2];
sx q[2];
rz(-1.0641629) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.74789819) q[1];
sx q[1];
rz(-1.6130578) q[1];
sx q[1];
rz(1.7903922) q[1];
rz(1.140959) q[3];
sx q[3];
rz(-1.4996734) q[3];
sx q[3];
rz(-1.5444559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0685136) q[2];
sx q[2];
rz(-1.2129236) q[2];
sx q[2];
rz(-0.53654137) q[2];
rz(2.1327298) q[3];
sx q[3];
rz(-0.77102414) q[3];
sx q[3];
rz(2.3276276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5040078) q[0];
sx q[0];
rz(-1.3115839) q[0];
sx q[0];
rz(0.017039321) q[0];
rz(-2.0783966) q[1];
sx q[1];
rz(-2.3737962) q[1];
sx q[1];
rz(-2.1868736) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5009618) q[0];
sx q[0];
rz(-2.9056773) q[0];
sx q[0];
rz(-2.9032034) q[0];
rz(-2.0644215) q[2];
sx q[2];
rz(-0.33702484) q[2];
sx q[2];
rz(-2.9309907) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8361289) q[1];
sx q[1];
rz(-1.3333798) q[1];
sx q[1];
rz(-1.8029455) q[1];
rz(-pi) q[2];
rz(3.122284) q[3];
sx q[3];
rz(-1.6240162) q[3];
sx q[3];
rz(-1.9198708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.22975989) q[2];
sx q[2];
rz(-0.91270295) q[2];
sx q[2];
rz(-2.0274053) q[2];
rz(-2.0071425) q[3];
sx q[3];
rz(-2.9454234) q[3];
sx q[3];
rz(-2.9451059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5512307) q[0];
sx q[0];
rz(-0.95142618) q[0];
sx q[0];
rz(2.9587342) q[0];
rz(3.0602449) q[1];
sx q[1];
rz(-2.3902939) q[1];
sx q[1];
rz(-0.61838165) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.189181) q[0];
sx q[0];
rz(-0.78972406) q[0];
sx q[0];
rz(3.038622) q[0];
x q[1];
rz(2.9419627) q[2];
sx q[2];
rz(-0.77130328) q[2];
sx q[2];
rz(-0.16801258) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.8733682) q[1];
sx q[1];
rz(-1.7098688) q[1];
sx q[1];
rz(1.2285608) q[1];
rz(-pi) q[2];
rz(-2.2719273) q[3];
sx q[3];
rz(-2.0732862) q[3];
sx q[3];
rz(-2.0521856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0028093) q[2];
sx q[2];
rz(-1.323779) q[2];
sx q[2];
rz(0.36661026) q[2];
rz(-1.9256516) q[3];
sx q[3];
rz(-1.9462908) q[3];
sx q[3];
rz(0.6634357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4902896) q[0];
sx q[0];
rz(-1.2325352) q[0];
sx q[0];
rz(2.41462) q[0];
rz(0.90826774) q[1];
sx q[1];
rz(-2.5636702) q[1];
sx q[1];
rz(1.04331) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7842877) q[0];
sx q[0];
rz(-2.4803964) q[0];
sx q[0];
rz(-3.1327598) q[0];
x q[1];
rz(1.7264446) q[2];
sx q[2];
rz(-2.7619432) q[2];
sx q[2];
rz(1.8774892) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9227495) q[1];
sx q[1];
rz(-2.6237539) q[1];
sx q[1];
rz(1.0475558) q[1];
rz(1.0087541) q[3];
sx q[3];
rz(-1.8971473) q[3];
sx q[3];
rz(-0.98451738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4116481) q[2];
sx q[2];
rz(-0.32117716) q[2];
sx q[2];
rz(-1.8357065) q[2];
rz(-0.11792396) q[3];
sx q[3];
rz(-0.4685466) q[3];
sx q[3];
rz(1.0606631) q[3];
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
rz(-2.1082728) q[0];
sx q[0];
rz(-0.98618996) q[0];
sx q[0];
rz(-1.5340075) q[0];
rz(-2.8458505) q[1];
sx q[1];
rz(-1.5402126) q[1];
sx q[1];
rz(1.6827513) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0223607) q[0];
sx q[0];
rz(-1.6408477) q[0];
sx q[0];
rz(-2.3168867) q[0];
rz(-pi) q[1];
rz(-0.88250156) q[2];
sx q[2];
rz(-2.6771328) q[2];
sx q[2];
rz(1.2683753) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.69469317) q[1];
sx q[1];
rz(-2.2647479) q[1];
sx q[1];
rz(0.25687053) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2269873) q[3];
sx q[3];
rz(-0.60887486) q[3];
sx q[3];
rz(-0.46907779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38833388) q[2];
sx q[2];
rz(-2.7497141) q[2];
sx q[2];
rz(-2.4957116) q[2];
rz(1.215747) q[3];
sx q[3];
rz(-1.682621) q[3];
sx q[3];
rz(1.7997883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.129824) q[0];
sx q[0];
rz(-2.3079066) q[0];
sx q[0];
rz(0.60761333) q[0];
rz(-1.9629078) q[1];
sx q[1];
rz(-1.8557529) q[1];
sx q[1];
rz(-0.36373055) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5967412) q[0];
sx q[0];
rz(-0.57692301) q[0];
sx q[0];
rz(2.0175009) q[0];
x q[1];
rz(-1.9035788) q[2];
sx q[2];
rz(-1.6146891) q[2];
sx q[2];
rz(1.8157168) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1235001) q[1];
sx q[1];
rz(-2.315633) q[1];
sx q[1];
rz(2.9781746) q[1];
rz(1.5341982) q[3];
sx q[3];
rz(-0.59662899) q[3];
sx q[3];
rz(-2.6109909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.58934775) q[2];
sx q[2];
rz(-0.10422464) q[2];
sx q[2];
rz(1.9488526) q[2];
rz(-2.4097811) q[3];
sx q[3];
rz(-1.7164427) q[3];
sx q[3];
rz(1.4881136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28855395) q[0];
sx q[0];
rz(-2.9712501) q[0];
sx q[0];
rz(-1.6188251) q[0];
rz(1.4672) q[1];
sx q[1];
rz(-1.3102945) q[1];
sx q[1];
rz(-2.4640962) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2427432) q[0];
sx q[0];
rz(-1.5206608) q[0];
sx q[0];
rz(2.2671347) q[0];
rz(-0.23394312) q[2];
sx q[2];
rz(-1.3877227) q[2];
sx q[2];
rz(-1.2884566) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4851393) q[1];
sx q[1];
rz(-2.5514126) q[1];
sx q[1];
rz(-0.82832576) q[1];
x q[2];
rz(-0.17685276) q[3];
sx q[3];
rz(-0.74027432) q[3];
sx q[3];
rz(0.76903421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.39685321) q[2];
sx q[2];
rz(-0.91529673) q[2];
sx q[2];
rz(-1.8358561) q[2];
rz(0.33852494) q[3];
sx q[3];
rz(-1.1208231) q[3];
sx q[3];
rz(2.4566076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5386388) q[0];
sx q[0];
rz(-0.21553497) q[0];
sx q[0];
rz(2.9576874) q[0];
rz(1.6869847) q[1];
sx q[1];
rz(-2.589476) q[1];
sx q[1];
rz(-2.6457381) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60431474) q[0];
sx q[0];
rz(-1.0838944) q[0];
sx q[0];
rz(-2.8394798) q[0];
rz(2.3073436) q[2];
sx q[2];
rz(-2.1733892) q[2];
sx q[2];
rz(1.955223) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.42577463) q[1];
sx q[1];
rz(-1.1905021) q[1];
sx q[1];
rz(-2.7575995) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3794231) q[3];
sx q[3];
rz(-0.56063214) q[3];
sx q[3];
rz(2.3440685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.054606525) q[2];
sx q[2];
rz(-0.84422529) q[2];
sx q[2];
rz(1.1620713) q[2];
rz(-1.453513) q[3];
sx q[3];
rz(-0.96265692) q[3];
sx q[3];
rz(-1.2801722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7206955) q[0];
sx q[0];
rz(-1.1524042) q[0];
sx q[0];
rz(-2.5757117) q[0];
rz(-0.3903009) q[1];
sx q[1];
rz(-1.5719599) q[1];
sx q[1];
rz(1.3077259) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5222675) q[0];
sx q[0];
rz(-0.59460708) q[0];
sx q[0];
rz(-0.99647605) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9487319) q[2];
sx q[2];
rz(-1.3728752) q[2];
sx q[2];
rz(1.6557616) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.81586397) q[1];
sx q[1];
rz(-0.96885175) q[1];
sx q[1];
rz(1.25156) q[1];
rz(-pi) q[2];
rz(-1.4212178) q[3];
sx q[3];
rz(-1.4523376) q[3];
sx q[3];
rz(1.3988023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2719443) q[2];
sx q[2];
rz(-2.7248236) q[2];
sx q[2];
rz(-2.1040253) q[2];
rz(1.8218254) q[3];
sx q[3];
rz(-2.3128553) q[3];
sx q[3];
rz(0.64798361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2607525) q[0];
sx q[0];
rz(-2.1639731) q[0];
sx q[0];
rz(0.64224893) q[0];
rz(1.3308659) q[1];
sx q[1];
rz(-0.86527491) q[1];
sx q[1];
rz(2.9790402) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7851622) q[0];
sx q[0];
rz(-0.9097865) q[0];
sx q[0];
rz(-0.50161538) q[0];
x q[1];
rz(-0.37940782) q[2];
sx q[2];
rz(-1.9072633) q[2];
sx q[2];
rz(3.0010146) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6813415) q[1];
sx q[1];
rz(-0.5149018) q[1];
sx q[1];
rz(1.9327362) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0216072) q[3];
sx q[3];
rz(-2.5549742) q[3];
sx q[3];
rz(-0.86931673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6315397) q[2];
sx q[2];
rz(-1.8659464) q[2];
sx q[2];
rz(2.1007382) q[2];
rz(0.96796525) q[3];
sx q[3];
rz(-2.0716045) q[3];
sx q[3];
rz(0.66106838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80617245) q[0];
sx q[0];
rz(-1.5763043) q[0];
sx q[0];
rz(3.0446654) q[0];
rz(-0.23282911) q[1];
sx q[1];
rz(-2.1131344) q[1];
sx q[1];
rz(0.0075385787) q[1];
rz(2.0816879) q[2];
sx q[2];
rz(-2.2507406) q[2];
sx q[2];
rz(3.0377664) q[2];
rz(2.6951058) q[3];
sx q[3];
rz(-0.75110186) q[3];
sx q[3];
rz(2.5121381) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
