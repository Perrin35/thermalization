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
rz(-0.023802726) q[0];
sx q[0];
rz(4.2476141) q[0];
sx q[0];
rz(10.201693) q[0];
rz(-1.7493526) q[1];
sx q[1];
rz(-1.8267781) q[1];
sx q[1];
rz(-2.1652752) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1084749) q[0];
sx q[0];
rz(-0.43370789) q[0];
sx q[0];
rz(-1.2483622) q[0];
rz(-pi) q[1];
rz(-1.5037276) q[2];
sx q[2];
rz(-2.133495) q[2];
sx q[2];
rz(-2.2864311) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9982517) q[1];
sx q[1];
rz(-0.51330459) q[1];
sx q[1];
rz(0.99350382) q[1];
x q[2];
rz(1.5209274) q[3];
sx q[3];
rz(-0.5736151) q[3];
sx q[3];
rz(-0.22465868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6355847) q[2];
sx q[2];
rz(-0.053746544) q[2];
sx q[2];
rz(-0.40679833) q[2];
rz(-2.9721416) q[3];
sx q[3];
rz(-0.5295161) q[3];
sx q[3];
rz(1.0725526) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2605543) q[0];
sx q[0];
rz(-0.23242234) q[0];
sx q[0];
rz(-3.1378003) q[0];
rz(-3.0637528) q[1];
sx q[1];
rz(-0.65996116) q[1];
sx q[1];
rz(-0.30581623) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90214848) q[0];
sx q[0];
rz(-1.3864707) q[0];
sx q[0];
rz(2.091616) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7051093) q[2];
sx q[2];
rz(-0.82342734) q[2];
sx q[2];
rz(-3.1253377) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7485436) q[1];
sx q[1];
rz(-0.83213388) q[1];
sx q[1];
rz(2.119675) q[1];
x q[2];
rz(3.0806957) q[3];
sx q[3];
rz(-1.6346667) q[3];
sx q[3];
rz(-2.8012432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99016142) q[2];
sx q[2];
rz(-1.3913245) q[2];
sx q[2];
rz(-2.4988417) q[2];
rz(-0.07240545) q[3];
sx q[3];
rz(-1.0591155) q[3];
sx q[3];
rz(-1.36093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0533503) q[0];
sx q[0];
rz(-2.998816) q[0];
sx q[0];
rz(2.9852168) q[0];
rz(3.1001672) q[1];
sx q[1];
rz(-2.5138469) q[1];
sx q[1];
rz(1.5511537) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0254733) q[0];
sx q[0];
rz(-2.1120694) q[0];
sx q[0];
rz(1.3295637) q[0];
x q[1];
rz(0.088038283) q[2];
sx q[2];
rz(-0.65795846) q[2];
sx q[2];
rz(-0.41706271) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.52122766) q[1];
sx q[1];
rz(-1.4583734) q[1];
sx q[1];
rz(-0.53037723) q[1];
x q[2];
rz(0.87167344) q[3];
sx q[3];
rz(-2.1101885) q[3];
sx q[3];
rz(1.9339069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7330043) q[2];
sx q[2];
rz(-0.53885794) q[2];
sx q[2];
rz(0.020922529) q[2];
rz(-2.9544592) q[3];
sx q[3];
rz(-0.20550607) q[3];
sx q[3];
rz(3.0278897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1784096) q[0];
sx q[0];
rz(-0.33247501) q[0];
sx q[0];
rz(-0.43854976) q[0];
rz(-1.6167538) q[1];
sx q[1];
rz(-0.33477819) q[1];
sx q[1];
rz(2.8964892) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.995718) q[0];
sx q[0];
rz(-0.85332131) q[0];
sx q[0];
rz(3.0307253) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39056492) q[2];
sx q[2];
rz(-0.3949983) q[2];
sx q[2];
rz(-2.4633138) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.41293834) q[1];
sx q[1];
rz(-1.0333583) q[1];
sx q[1];
rz(2.9793903) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3880344) q[3];
sx q[3];
rz(-1.8055969) q[3];
sx q[3];
rz(-0.59111349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.54030067) q[2];
sx q[2];
rz(-0.43593323) q[2];
sx q[2];
rz(0.33176804) q[2];
rz(0.48745421) q[3];
sx q[3];
rz(-2.09477) q[3];
sx q[3];
rz(-0.99307466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29193923) q[0];
sx q[0];
rz(-1.6878457) q[0];
sx q[0];
rz(0.77350235) q[0];
rz(-1.1812814) q[1];
sx q[1];
rz(-0.14207323) q[1];
sx q[1];
rz(-1.7519417) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4163602) q[0];
sx q[0];
rz(-1.9342039) q[0];
sx q[0];
rz(-2.711722) q[0];
rz(-pi) q[1];
rz(-0.85948617) q[2];
sx q[2];
rz(-2.0469249) q[2];
sx q[2];
rz(1.7993594) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.70848318) q[1];
sx q[1];
rz(-2.6603087) q[1];
sx q[1];
rz(2.9631056) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64719154) q[3];
sx q[3];
rz(-2.8464937) q[3];
sx q[3];
rz(-1.6012675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2191849) q[2];
sx q[2];
rz(-1.9318523) q[2];
sx q[2];
rz(-0.47214559) q[2];
rz(-1.8426497) q[3];
sx q[3];
rz(-1.2932212) q[3];
sx q[3];
rz(-0.81056547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2209114) q[0];
sx q[0];
rz(-0.26316106) q[0];
sx q[0];
rz(-2.8780908) q[0];
rz(1.1031411) q[1];
sx q[1];
rz(-1.8270854) q[1];
sx q[1];
rz(2.7679494) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3220396) q[0];
sx q[0];
rz(-3.0470938) q[0];
sx q[0];
rz(-2.2631133) q[0];
x q[1];
rz(0.12219001) q[2];
sx q[2];
rz(-0.8475248) q[2];
sx q[2];
rz(3.0360589) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9615104) q[1];
sx q[1];
rz(-2.1282548) q[1];
sx q[1];
rz(0.37661676) q[1];
rz(-pi) q[2];
rz(0.5881891) q[3];
sx q[3];
rz(-1.08687) q[3];
sx q[3];
rz(1.5744792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0282447) q[2];
sx q[2];
rz(-2.9714606) q[2];
sx q[2];
rz(0.55383468) q[2];
rz(1.3977741) q[3];
sx q[3];
rz(-2.5388986) q[3];
sx q[3];
rz(-2.8288614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5572307) q[0];
sx q[0];
rz(-2.1350242) q[0];
sx q[0];
rz(-1.2297909) q[0];
rz(2.9025485) q[1];
sx q[1];
rz(-1.6300647) q[1];
sx q[1];
rz(0.30034932) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8431138) q[0];
sx q[0];
rz(-1.686578) q[0];
sx q[0];
rz(-2.9779469) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0763758) q[2];
sx q[2];
rz(-0.9305312) q[2];
sx q[2];
rz(2.6312021) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9253941) q[1];
sx q[1];
rz(-1.5860737) q[1];
sx q[1];
rz(3.141204) q[1];
rz(-pi) q[2];
x q[2];
rz(0.024552931) q[3];
sx q[3];
rz(-0.85806393) q[3];
sx q[3];
rz(-0.36125444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0099237) q[2];
sx q[2];
rz(-1.6841623) q[2];
sx q[2];
rz(0.24492502) q[2];
rz(-0.51982546) q[3];
sx q[3];
rz(-0.86331415) q[3];
sx q[3];
rz(-0.68827099) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7898665) q[0];
sx q[0];
rz(-2.7930197) q[0];
sx q[0];
rz(-1.9785731) q[0];
rz(-3.0746958) q[1];
sx q[1];
rz(-1.6480548) q[1];
sx q[1];
rz(-1.012872) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.430535) q[0];
sx q[0];
rz(-1.0037046) q[0];
sx q[0];
rz(2.5398769) q[0];
rz(1.0032907) q[2];
sx q[2];
rz(-1.7219647) q[2];
sx q[2];
rz(1.0790107) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7808395) q[1];
sx q[1];
rz(-1.2560802) q[1];
sx q[1];
rz(0.34784045) q[1];
rz(1.1539641) q[3];
sx q[3];
rz(-0.92769054) q[3];
sx q[3];
rz(-2.7288306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0248727) q[2];
sx q[2];
rz(-2.1634384) q[2];
sx q[2];
rz(-0.32279521) q[2];
rz(0.60574496) q[3];
sx q[3];
rz(-2.3492458) q[3];
sx q[3];
rz(-0.34887031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.087273) q[0];
sx q[0];
rz(-0.1467341) q[0];
sx q[0];
rz(0.012454575) q[0];
rz(2.3948578) q[1];
sx q[1];
rz(-2.2181999) q[1];
sx q[1];
rz(-2.8616203) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075172193) q[0];
sx q[0];
rz(-3.018258) q[0];
sx q[0];
rz(-1.146011) q[0];
rz(-1.1367646) q[2];
sx q[2];
rz(-1.869259) q[2];
sx q[2];
rz(0.091574319) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.88290434) q[1];
sx q[1];
rz(-1.2637485) q[1];
sx q[1];
rz(-2.3767002) q[1];
rz(-pi) q[2];
rz(-2.1567508) q[3];
sx q[3];
rz(-2.7142314) q[3];
sx q[3];
rz(-2.7387184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3296457) q[2];
sx q[2];
rz(-0.75355607) q[2];
sx q[2];
rz(2.9296056) q[2];
rz(0.82344615) q[3];
sx q[3];
rz(-1.6971089) q[3];
sx q[3];
rz(0.25920355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9987746) q[0];
sx q[0];
rz(-0.057567216) q[0];
sx q[0];
rz(-2.4488191) q[0];
rz(-0.57299262) q[1];
sx q[1];
rz(-1.3447821) q[1];
sx q[1];
rz(-2.7105892) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13700542) q[0];
sx q[0];
rz(-2.1218897) q[0];
sx q[0];
rz(-1.5020834) q[0];
rz(-pi) q[1];
rz(0.60811483) q[2];
sx q[2];
rz(-0.69497847) q[2];
sx q[2];
rz(2.8335477) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0242053) q[1];
sx q[1];
rz(-2.3820602) q[1];
sx q[1];
rz(2.7717436) q[1];
rz(-pi) q[2];
rz(-1.0968571) q[3];
sx q[3];
rz(-2.5383679) q[3];
sx q[3];
rz(-2.4639377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7196322) q[2];
sx q[2];
rz(-2.8675291) q[2];
sx q[2];
rz(-0.56023041) q[2];
rz(-2.6719921) q[3];
sx q[3];
rz(-2.738939) q[3];
sx q[3];
rz(0.71389055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8568759) q[0];
sx q[0];
rz(-1.7354043) q[0];
sx q[0];
rz(2.0866557) q[0];
rz(-0.84125413) q[1];
sx q[1];
rz(-2.0396736) q[1];
sx q[1];
rz(3.094818) q[1];
rz(-1.9404491) q[2];
sx q[2];
rz(-1.2751725) q[2];
sx q[2];
rz(-1.0775492) q[2];
rz(1.6056521) q[3];
sx q[3];
rz(-0.65924725) q[3];
sx q[3];
rz(1.02871) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
