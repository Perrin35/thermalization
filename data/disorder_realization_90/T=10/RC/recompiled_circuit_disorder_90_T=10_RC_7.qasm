OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7251627) q[0];
sx q[0];
rz(-3.0017612) q[0];
sx q[0];
rz(-0.60959417) q[0];
rz(-5.6145515) q[1];
sx q[1];
rz(0.86548391) q[1];
sx q[1];
rz(15.794985) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60458175) q[0];
sx q[0];
rz(-0.7365948) q[0];
sx q[0];
rz(2.1405311) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1331698) q[2];
sx q[2];
rz(-0.77169092) q[2];
sx q[2];
rz(2.0193677) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0793003) q[1];
sx q[1];
rz(-1.0907409) q[1];
sx q[1];
rz(0.15602195) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.617241) q[3];
sx q[3];
rz(-1.9825476) q[3];
sx q[3];
rz(2.9573033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0027851) q[2];
sx q[2];
rz(-1.3802718) q[2];
sx q[2];
rz(-0.37386093) q[2];
rz(2.8047681) q[3];
sx q[3];
rz(-1.5461494) q[3];
sx q[3];
rz(0.22836223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(1.2574629) q[0];
sx q[0];
rz(-0.3733491) q[0];
sx q[0];
rz(1.194838) q[0];
rz(-3.0589814) q[1];
sx q[1];
rz(-1.9742842) q[1];
sx q[1];
rz(-3.1412178) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9449687) q[0];
sx q[0];
rz(-1.9163016) q[0];
sx q[0];
rz(1.1254805) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20034321) q[2];
sx q[2];
rz(-1.3140972) q[2];
sx q[2];
rz(0.30470195) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5217168) q[1];
sx q[1];
rz(-2.4881425) q[1];
sx q[1];
rz(0.39342777) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5156636) q[3];
sx q[3];
rz(-1.1109567) q[3];
sx q[3];
rz(2.9210747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.80883819) q[2];
sx q[2];
rz(-0.19503441) q[2];
sx q[2];
rz(0.17671281) q[2];
rz(0.79408944) q[3];
sx q[3];
rz(-0.77787557) q[3];
sx q[3];
rz(2.5879522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90213838) q[0];
sx q[0];
rz(-0.70394009) q[0];
sx q[0];
rz(-0.40476558) q[0];
rz(-1.2813214) q[1];
sx q[1];
rz(-0.57360137) q[1];
sx q[1];
rz(1.8331029) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.20298) q[0];
sx q[0];
rz(-2.4327607) q[0];
sx q[0];
rz(-2.3908486) q[0];
x q[1];
rz(-0.53738014) q[2];
sx q[2];
rz(-1.3668622) q[2];
sx q[2];
rz(0.98762074) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.96924671) q[1];
sx q[1];
rz(-1.8241276) q[1];
sx q[1];
rz(0.84448703) q[1];
rz(-pi) q[2];
rz(-0.43867302) q[3];
sx q[3];
rz(-1.6524501) q[3];
sx q[3];
rz(-0.7338394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.219316) q[2];
sx q[2];
rz(-0.52307659) q[2];
sx q[2];
rz(-3.1211839) q[2];
rz(2.0698047) q[3];
sx q[3];
rz(-2.0139549) q[3];
sx q[3];
rz(2.4782457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9999009) q[0];
sx q[0];
rz(-0.90943709) q[0];
sx q[0];
rz(0.91598696) q[0];
rz(-2.6782716) q[1];
sx q[1];
rz(-2.0756192) q[1];
sx q[1];
rz(-1.0571009) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6592641) q[0];
sx q[0];
rz(-1.5963012) q[0];
sx q[0];
rz(0.43735023) q[0];
rz(-2.1787203) q[2];
sx q[2];
rz(-2.0553556) q[2];
sx q[2];
rz(1.0819266) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5406487) q[1];
sx q[1];
rz(-1.800866) q[1];
sx q[1];
rz(0.01407108) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3841342) q[3];
sx q[3];
rz(-1.7372903) q[3];
sx q[3];
rz(-1.43653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2477734) q[2];
sx q[2];
rz(-0.63988581) q[2];
sx q[2];
rz(-2.7988953) q[2];
rz(1.6977067) q[3];
sx q[3];
rz(-0.87564898) q[3];
sx q[3];
rz(0.7152344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3559568) q[0];
sx q[0];
rz(-2.5781093) q[0];
sx q[0];
rz(0.086439565) q[0];
rz(-1.3899639) q[1];
sx q[1];
rz(-0.93170634) q[1];
sx q[1];
rz(0.46868971) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5422573) q[0];
sx q[0];
rz(-2.242803) q[0];
sx q[0];
rz(2.6319648) q[0];
rz(-pi) q[1];
rz(2.8332963) q[2];
sx q[2];
rz(-1.397965) q[2];
sx q[2];
rz(2.1413213) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9710755) q[1];
sx q[1];
rz(-1.3689539) q[1];
sx q[1];
rz(1.7032196) q[1];
rz(-1.4851941) q[3];
sx q[3];
rz(-2.3833131) q[3];
sx q[3];
rz(-3.0183834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.1420574) q[2];
sx q[2];
rz(-2.2613328) q[2];
sx q[2];
rz(-2.2772677) q[2];
rz(0.21720973) q[3];
sx q[3];
rz(-0.61621284) q[3];
sx q[3];
rz(-0.29278452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42751673) q[0];
sx q[0];
rz(-1.7345411) q[0];
sx q[0];
rz(-0.19700225) q[0];
rz(-1.7794094) q[1];
sx q[1];
rz(-1.660659) q[1];
sx q[1];
rz(-0.33624712) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7173548) q[0];
sx q[0];
rz(-2.3399118) q[0];
sx q[0];
rz(-2.5108811) q[0];
rz(0.98724987) q[2];
sx q[2];
rz(-1.0580214) q[2];
sx q[2];
rz(3.0836881) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.044144883) q[1];
sx q[1];
rz(-2.1598585) q[1];
sx q[1];
rz(1.6277905) q[1];
rz(-pi) q[2];
rz(-1.398596) q[3];
sx q[3];
rz(-2.2599054) q[3];
sx q[3];
rz(0.25758753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6062935) q[2];
sx q[2];
rz(-1.4761816) q[2];
sx q[2];
rz(0.84632787) q[2];
rz(1.2396631) q[3];
sx q[3];
rz(-1.2318434) q[3];
sx q[3];
rz(-3.1404176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7476615) q[0];
sx q[0];
rz(-1.0378391) q[0];
sx q[0];
rz(2.8826662) q[0];
rz(1.3461643) q[1];
sx q[1];
rz(-1.3849473) q[1];
sx q[1];
rz(1.0940201) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0712229) q[0];
sx q[0];
rz(-0.56814146) q[0];
sx q[0];
rz(2.2413261) q[0];
rz(1.3482434) q[2];
sx q[2];
rz(-1.5587274) q[2];
sx q[2];
rz(-2.0633069) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.5309696) q[1];
sx q[1];
rz(-1.5986643) q[1];
sx q[1];
rz(-1.665297) q[1];
rz(1.5488946) q[3];
sx q[3];
rz(-1.4727482) q[3];
sx q[3];
rz(-2.736562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9672433) q[2];
sx q[2];
rz(-1.7847585) q[2];
sx q[2];
rz(0.32067498) q[2];
rz(-0.43618068) q[3];
sx q[3];
rz(-0.47302055) q[3];
sx q[3];
rz(0.44803739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2748579) q[0];
sx q[0];
rz(-0.47645706) q[0];
sx q[0];
rz(2.136769) q[0];
rz(2.5514305) q[1];
sx q[1];
rz(-0.90548038) q[1];
sx q[1];
rz(0.80387962) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.852254) q[0];
sx q[0];
rz(-2.5544871) q[0];
sx q[0];
rz(-0.071285204) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5875823) q[2];
sx q[2];
rz(-1.4129352) q[2];
sx q[2];
rz(0.49288921) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.67919532) q[1];
sx q[1];
rz(-0.089226626) q[1];
sx q[1];
rz(-1.4560844) q[1];
rz(-1.6611093) q[3];
sx q[3];
rz(-0.90925018) q[3];
sx q[3];
rz(-2.6722398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8089495) q[2];
sx q[2];
rz(-0.84471622) q[2];
sx q[2];
rz(-1.0104898) q[2];
rz(0.60339749) q[3];
sx q[3];
rz(-1.6167275) q[3];
sx q[3];
rz(-2.5914014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5378961) q[0];
sx q[0];
rz(-1.2768856) q[0];
sx q[0];
rz(-2.7340775) q[0];
rz(-0.28911668) q[1];
sx q[1];
rz(-2.018785) q[1];
sx q[1];
rz(-0.75072748) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4780316) q[0];
sx q[0];
rz(-2.7285517) q[0];
sx q[0];
rz(2.3703299) q[0];
rz(-pi) q[1];
rz(1.6523916) q[2];
sx q[2];
rz(-1.8494693) q[2];
sx q[2];
rz(0.96688731) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0093065) q[1];
sx q[1];
rz(-2.7704151) q[1];
sx q[1];
rz(2.1585141) q[1];
x q[2];
rz(1.5446072) q[3];
sx q[3];
rz(-1.8720172) q[3];
sx q[3];
rz(-1.6645886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0728545) q[2];
sx q[2];
rz(-0.72454238) q[2];
sx q[2];
rz(-1.9753974) q[2];
rz(-1.7769622) q[3];
sx q[3];
rz(-0.77562538) q[3];
sx q[3];
rz(-3.1197746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5584548) q[0];
sx q[0];
rz(-2.3151509) q[0];
sx q[0];
rz(1.3903842) q[0];
rz(-2.8109) q[1];
sx q[1];
rz(-2.3762517) q[1];
sx q[1];
rz(-1.6814544) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18171039) q[0];
sx q[0];
rz(-1.0707756) q[0];
sx q[0];
rz(2.659003) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5286469) q[2];
sx q[2];
rz(-0.50439207) q[2];
sx q[2];
rz(1.1500037) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9863661) q[1];
sx q[1];
rz(-0.44253293) q[1];
sx q[1];
rz(0.44154422) q[1];
x q[2];
rz(-0.86942418) q[3];
sx q[3];
rz(-0.93360177) q[3];
sx q[3];
rz(0.77995342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.34974393) q[2];
sx q[2];
rz(-1.8169553) q[2];
sx q[2];
rz(-1.6798518) q[2];
rz(2.0215624) q[3];
sx q[3];
rz(-2.5199065) q[3];
sx q[3];
rz(0.45564836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5158952) q[0];
sx q[0];
rz(-1.6202171) q[0];
sx q[0];
rz(1.1465999) q[0];
rz(1.3810146) q[1];
sx q[1];
rz(-1.2881423) q[1];
sx q[1];
rz(-1.2013411) q[1];
rz(-0.78787055) q[2];
sx q[2];
rz(-1.315016) q[2];
sx q[2];
rz(-1.6223326) q[2];
rz(-2.667726) q[3];
sx q[3];
rz(-2.4110473) q[3];
sx q[3];
rz(-1.9163781) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
