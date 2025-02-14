OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5918936) q[0];
sx q[0];
rz(-2.0285719) q[0];
sx q[0];
rz(-0.23166238) q[0];
rz(0.35056937) q[1];
sx q[1];
rz(6.0729519) q[1];
sx q[1];
rz(8.6624866) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9755154) q[0];
sx q[0];
rz(-0.78801256) q[0];
sx q[0];
rz(-0.60933526) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1074802) q[2];
sx q[2];
rz(-1.1321403) q[2];
sx q[2];
rz(0.92213501) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8787391) q[1];
sx q[1];
rz(-1.384842) q[1];
sx q[1];
rz(-1.6270217) q[1];
x q[2];
rz(-2.2985994) q[3];
sx q[3];
rz(-0.61474934) q[3];
sx q[3];
rz(0.45045567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1594557) q[2];
sx q[2];
rz(-0.60672131) q[2];
sx q[2];
rz(-0.43178001) q[2];
rz(-0.86709705) q[3];
sx q[3];
rz(-0.95621395) q[3];
sx q[3];
rz(1.506327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1113488) q[0];
sx q[0];
rz(-2.3115277) q[0];
sx q[0];
rz(1.9507971) q[0];
rz(-1.8973154) q[1];
sx q[1];
rz(-2.0121274) q[1];
sx q[1];
rz(0.58070374) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8180926) q[0];
sx q[0];
rz(-0.43002263) q[0];
sx q[0];
rz(-1.9097206) q[0];
rz(-2.2207308) q[2];
sx q[2];
rz(-2.5061786) q[2];
sx q[2];
rz(0.17104761) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8088219) q[1];
sx q[1];
rz(-0.94478196) q[1];
sx q[1];
rz(1.9653734) q[1];
rz(-pi) q[2];
rz(-1.034581) q[3];
sx q[3];
rz(-2.381122) q[3];
sx q[3];
rz(-1.5813259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1397436) q[2];
sx q[2];
rz(-1.7364343) q[2];
sx q[2];
rz(-1.5025567) q[2];
rz(2.9338845) q[3];
sx q[3];
rz(-1.3108871) q[3];
sx q[3];
rz(-2.1343855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.945632) q[0];
sx q[0];
rz(-1.0769083) q[0];
sx q[0];
rz(2.9929602) q[0];
rz(2.0678988) q[1];
sx q[1];
rz(-0.37207347) q[1];
sx q[1];
rz(-0.78280848) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4769094) q[0];
sx q[0];
rz(-2.9626155) q[0];
sx q[0];
rz(-1.2577235) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6086968) q[2];
sx q[2];
rz(-1.6072465) q[2];
sx q[2];
rz(-0.25668555) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9298934) q[1];
sx q[1];
rz(-0.5553588) q[1];
sx q[1];
rz(-2.7564605) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25198307) q[3];
sx q[3];
rz(-2.3445233) q[3];
sx q[3];
rz(-1.0154203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.17915501) q[2];
sx q[2];
rz(-2.16733) q[2];
sx q[2];
rz(1.1243189) q[2];
rz(3.1149241) q[3];
sx q[3];
rz(-2.4494438) q[3];
sx q[3];
rz(2.8587225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35793316) q[0];
sx q[0];
rz(-1.986035) q[0];
sx q[0];
rz(1.5643157) q[0];
rz(-1.0464926) q[1];
sx q[1];
rz(-1.0341897) q[1];
sx q[1];
rz(2.0139205) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8148147) q[0];
sx q[0];
rz(-2.5501857) q[0];
sx q[0];
rz(3.1230984) q[0];
x q[1];
rz(1.7972184) q[2];
sx q[2];
rz(-0.34219301) q[2];
sx q[2];
rz(0.36384091) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8279523) q[1];
sx q[1];
rz(-1.2544187) q[1];
sx q[1];
rz(-0.72319855) q[1];
rz(1.4521056) q[3];
sx q[3];
rz(-1.4726536) q[3];
sx q[3];
rz(-0.019358033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.54101888) q[2];
sx q[2];
rz(-1.8709196) q[2];
sx q[2];
rz(-0.18826558) q[2];
rz(-2.6514734) q[3];
sx q[3];
rz(-1.6907588) q[3];
sx q[3];
rz(1.8671487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1566496) q[0];
sx q[0];
rz(-1.5771414) q[0];
sx q[0];
rz(-1.6749325) q[0];
rz(-2.6690392) q[1];
sx q[1];
rz(-0.87635374) q[1];
sx q[1];
rz(2.1579425) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88681038) q[0];
sx q[0];
rz(-1.3894102) q[0];
sx q[0];
rz(-0.43116335) q[0];
rz(-2.7698414) q[2];
sx q[2];
rz(-2.6494827) q[2];
sx q[2];
rz(0.12002698) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77351219) q[1];
sx q[1];
rz(-1.0607914) q[1];
sx q[1];
rz(0.065392134) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9463169) q[3];
sx q[3];
rz(-0.28464475) q[3];
sx q[3];
rz(-1.1840905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5094362) q[2];
sx q[2];
rz(-1.7551273) q[2];
sx q[2];
rz(-2.7344088) q[2];
rz(-1.8699649) q[3];
sx q[3];
rz(-0.41912246) q[3];
sx q[3];
rz(2.4841888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7206409) q[0];
sx q[0];
rz(-1.9176418) q[0];
sx q[0];
rz(1.2366914) q[0];
rz(1.1300794) q[1];
sx q[1];
rz(-1.7944444) q[1];
sx q[1];
rz(2.0225661) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84270714) q[0];
sx q[0];
rz(-1.7241553) q[0];
sx q[0];
rz(-1.6394079) q[0];
rz(0.58216146) q[2];
sx q[2];
rz(-0.88554875) q[2];
sx q[2];
rz(0.2174046) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5608852) q[1];
sx q[1];
rz(-1.4262154) q[1];
sx q[1];
rz(-2.7156262) q[1];
rz(-pi) q[2];
rz(2.7009526) q[3];
sx q[3];
rz(-1.2088606) q[3];
sx q[3];
rz(-2.831865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.59469026) q[2];
sx q[2];
rz(-1.12744) q[2];
sx q[2];
rz(0.31000578) q[2];
rz(-1.7075432) q[3];
sx q[3];
rz(-1.6238345) q[3];
sx q[3];
rz(1.2448509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2731648) q[0];
sx q[0];
rz(-2.7719066) q[0];
sx q[0];
rz(-1.0721068) q[0];
rz(1.781069) q[1];
sx q[1];
rz(-0.37630263) q[1];
sx q[1];
rz(1.2881813) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38561197) q[0];
sx q[0];
rz(-0.98273425) q[0];
sx q[0];
rz(-0.62784348) q[0];
x q[1];
rz(2.6990876) q[2];
sx q[2];
rz(-1.0396084) q[2];
sx q[2];
rz(-2.1489904) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5406193) q[1];
sx q[1];
rz(-0.29635591) q[1];
sx q[1];
rz(-2.5286416) q[1];
x q[2];
rz(-0.76948072) q[3];
sx q[3];
rz(-1.2043287) q[3];
sx q[3];
rz(2.0539157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0897022) q[2];
sx q[2];
rz(-1.6334198) q[2];
sx q[2];
rz(1.7066329) q[2];
rz(-2.6658304) q[3];
sx q[3];
rz(-2.4475554) q[3];
sx q[3];
rz(-2.7580875) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7504904) q[0];
sx q[0];
rz(-0.35677156) q[0];
sx q[0];
rz(-1.179689) q[0];
rz(0.36243311) q[1];
sx q[1];
rz(-2.0795627) q[1];
sx q[1];
rz(-1.5865883) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3488551) q[0];
sx q[0];
rz(-2.7842583) q[0];
sx q[0];
rz(2.6941195) q[0];
x q[1];
rz(-2.4612975) q[2];
sx q[2];
rz(-1.9501424) q[2];
sx q[2];
rz(-0.85000402) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5240979) q[1];
sx q[1];
rz(-0.84071181) q[1];
sx q[1];
rz(2.1725656) q[1];
rz(-pi) q[2];
rz(-0.15896564) q[3];
sx q[3];
rz(-2.1329125) q[3];
sx q[3];
rz(-1.5288439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.72593752) q[2];
sx q[2];
rz(-1.1154117) q[2];
sx q[2];
rz(-1.1768781) q[2];
rz(-0.40500179) q[3];
sx q[3];
rz(-0.98086762) q[3];
sx q[3];
rz(0.39308959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28615752) q[0];
sx q[0];
rz(-3.0431008) q[0];
sx q[0];
rz(-1.806102) q[0];
rz(-2.2291741) q[1];
sx q[1];
rz(-1.4211979) q[1];
sx q[1];
rz(0.040249912) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9837058) q[0];
sx q[0];
rz(-1.5327546) q[0];
sx q[0];
rz(1.4020919) q[0];
rz(0.59800541) q[2];
sx q[2];
rz(-2.0705418) q[2];
sx q[2];
rz(2.6464484) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9973491) q[1];
sx q[1];
rz(-1.9895619) q[1];
sx q[1];
rz(2.1476157) q[1];
rz(-pi) q[2];
rz(-2.9222708) q[3];
sx q[3];
rz(-0.78028934) q[3];
sx q[3];
rz(-1.8530451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5361629) q[2];
sx q[2];
rz(-1.9998113) q[2];
sx q[2];
rz(-2.6119168) q[2];
rz(1.4179432) q[3];
sx q[3];
rz(-0.46897408) q[3];
sx q[3];
rz(-2.6403707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0025075992) q[0];
sx q[0];
rz(-1.7739828) q[0];
sx q[0];
rz(0.031877192) q[0];
rz(1.9335951) q[1];
sx q[1];
rz(-2.5937158) q[1];
sx q[1];
rz(0.30778232) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9986711) q[0];
sx q[0];
rz(-1.0772155) q[0];
sx q[0];
rz(-2.837985) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8372907) q[2];
sx q[2];
rz(-1.4478168) q[2];
sx q[2];
rz(2.3044555) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.59474241) q[1];
sx q[1];
rz(-1.062927) q[1];
sx q[1];
rz(2.2617729) q[1];
x q[2];
rz(2.4069837) q[3];
sx q[3];
rz(-0.19140581) q[3];
sx q[3];
rz(-3.1212487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6098183) q[2];
sx q[2];
rz(-2.4419624) q[2];
sx q[2];
rz(0.043665234) q[2];
rz(-1.4788491) q[3];
sx q[3];
rz(-0.50373977) q[3];
sx q[3];
rz(1.7536722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26535784) q[0];
sx q[0];
rz(-1.6725412) q[0];
sx q[0];
rz(-1.5681736) q[0];
rz(-3.1126032) q[1];
sx q[1];
rz(-2.034076) q[1];
sx q[1];
rz(-2.1709002) q[1];
rz(-0.12262298) q[2];
sx q[2];
rz(-1.7896252) q[2];
sx q[2];
rz(0.26364506) q[2];
rz(0.51359691) q[3];
sx q[3];
rz(-0.83181341) q[3];
sx q[3];
rz(0.60461525) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
