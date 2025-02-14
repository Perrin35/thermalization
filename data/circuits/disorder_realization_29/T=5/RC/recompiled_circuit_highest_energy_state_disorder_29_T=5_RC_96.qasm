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
rz(-0.40773243) q[0];
sx q[0];
rz(-1.9324349) q[0];
sx q[0];
rz(1.6785167) q[0];
rz(-2.8830124) q[1];
sx q[1];
rz(-1.2376031) q[1];
sx q[1];
rz(-3.0078476) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.84496) q[0];
sx q[0];
rz(-1.4136853) q[0];
sx q[0];
rz(-0.0060307302) q[0];
x q[1];
rz(2.1315246) q[2];
sx q[2];
rz(-1.0655996) q[2];
sx q[2];
rz(1.2743735) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7642149) q[1];
sx q[1];
rz(-2.1020319) q[1];
sx q[1];
rz(-1.9218535) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5695581) q[3];
sx q[3];
rz(-0.15286013) q[3];
sx q[3];
rz(-2.3521956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3818843) q[2];
sx q[2];
rz(-0.81177652) q[2];
sx q[2];
rz(1.3105357) q[2];
rz(-1.367502) q[3];
sx q[3];
rz(-2.01367) q[3];
sx q[3];
rz(3.1254356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97717184) q[0];
sx q[0];
rz(-1.9855969) q[0];
sx q[0];
rz(1.0154065) q[0];
rz(-2.7482391) q[1];
sx q[1];
rz(-1.669408) q[1];
sx q[1];
rz(-0.70158395) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4751101) q[0];
sx q[0];
rz(-0.85832867) q[0];
sx q[0];
rz(-0.11157596) q[0];
rz(1.6051454) q[2];
sx q[2];
rz(-0.30314244) q[2];
sx q[2];
rz(-0.042138635) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9394221) q[1];
sx q[1];
rz(-1.8864453) q[1];
sx q[1];
rz(3.062518) q[1];
rz(-pi) q[2];
rz(0.75357985) q[3];
sx q[3];
rz(-0.9155897) q[3];
sx q[3];
rz(1.1666949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1141438) q[2];
sx q[2];
rz(-1.0237209) q[2];
sx q[2];
rz(-2.0360428) q[2];
rz(-0.28087273) q[3];
sx q[3];
rz(-0.61087817) q[3];
sx q[3];
rz(-0.15225473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1059145) q[0];
sx q[0];
rz(-2.4690101) q[0];
sx q[0];
rz(1.0726844) q[0];
rz(0.013280344) q[1];
sx q[1];
rz(-1.599879) q[1];
sx q[1];
rz(-0.57659155) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7130797) q[0];
sx q[0];
rz(-1.858485) q[0];
sx q[0];
rz(-1.5031394) q[0];
rz(1.4378087) q[2];
sx q[2];
rz(-1.6536568) q[2];
sx q[2];
rz(-1.8216009) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.62200786) q[1];
sx q[1];
rz(-1.7652006) q[1];
sx q[1];
rz(2.0561051) q[1];
x q[2];
rz(-2.8956293) q[3];
sx q[3];
rz(-1.0268164) q[3];
sx q[3];
rz(1.1472052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.85792929) q[2];
sx q[2];
rz(-0.89202213) q[2];
sx q[2];
rz(2.499495) q[2];
rz(-0.56504956) q[3];
sx q[3];
rz(-0.27500209) q[3];
sx q[3];
rz(-0.16998418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2791486) q[0];
sx q[0];
rz(-1.2553517) q[0];
sx q[0];
rz(-0.24236648) q[0];
rz(-2.6077121) q[1];
sx q[1];
rz(-0.68478525) q[1];
sx q[1];
rz(-1.4885611) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80099312) q[0];
sx q[0];
rz(-1.6417608) q[0];
sx q[0];
rz(-2.5453955) q[0];
rz(-pi) q[1];
rz(-2.7127391) q[2];
sx q[2];
rz(-1.8629214) q[2];
sx q[2];
rz(-3.0512864) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6935255) q[1];
sx q[1];
rz(-0.84431813) q[1];
sx q[1];
rz(1.179715) q[1];
rz(-pi) q[2];
rz(-1.5307934) q[3];
sx q[3];
rz(-1.8703331) q[3];
sx q[3];
rz(-1.2114204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7786467) q[2];
sx q[2];
rz(-1.2641509) q[2];
sx q[2];
rz(-1.470559) q[2];
rz(0.10739022) q[3];
sx q[3];
rz(-1.8348179) q[3];
sx q[3];
rz(1.2764527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11112467) q[0];
sx q[0];
rz(-2.680439) q[0];
sx q[0];
rz(-1.9386468) q[0];
rz(-0.89490926) q[1];
sx q[1];
rz(-1.9590961) q[1];
sx q[1];
rz(-2.9295909) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1231736) q[0];
sx q[0];
rz(-0.97476649) q[0];
sx q[0];
rz(0.53431781) q[0];
x q[1];
rz(2.6190125) q[2];
sx q[2];
rz(-1.3088639) q[2];
sx q[2];
rz(-0.22522989) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.743269) q[1];
sx q[1];
rz(-1.124427) q[1];
sx q[1];
rz(-0.94029398) q[1];
rz(-pi) q[2];
rz(-1.5863442) q[3];
sx q[3];
rz(-0.55335303) q[3];
sx q[3];
rz(1.4912332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2682858) q[2];
sx q[2];
rz(-0.77750677) q[2];
sx q[2];
rz(3.0262465) q[2];
rz(2.1558732) q[3];
sx q[3];
rz(-1.7306354) q[3];
sx q[3];
rz(-0.43578291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2301521) q[0];
sx q[0];
rz(-1.6807012) q[0];
sx q[0];
rz(2.4440515) q[0];
rz(2.7279834) q[1];
sx q[1];
rz(-1.6802843) q[1];
sx q[1];
rz(0.70346171) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089982219) q[0];
sx q[0];
rz(-1.713366) q[0];
sx q[0];
rz(0.27326126) q[0];
rz(-pi) q[1];
x q[1];
rz(1.037961) q[2];
sx q[2];
rz(-2.46358) q[2];
sx q[2];
rz(0.17282669) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5825962) q[1];
sx q[1];
rz(-1.1038053) q[1];
sx q[1];
rz(-1.9860616) q[1];
rz(-pi) q[2];
x q[2];
rz(2.260842) q[3];
sx q[3];
rz(-0.38057571) q[3];
sx q[3];
rz(-2.5799023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.10282639) q[2];
sx q[2];
rz(-0.26831728) q[2];
sx q[2];
rz(1.212567) q[2];
rz(0.28193998) q[3];
sx q[3];
rz(-1.599267) q[3];
sx q[3];
rz(-0.50659242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7956227) q[0];
sx q[0];
rz(-2.3891734) q[0];
sx q[0];
rz(0.84683007) q[0];
rz(0.94547358) q[1];
sx q[1];
rz(-1.5739601) q[1];
sx q[1];
rz(-0.6792773) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69791171) q[0];
sx q[0];
rz(-2.6583932) q[0];
sx q[0];
rz(-2.9761821) q[0];
x q[1];
rz(2.5176454) q[2];
sx q[2];
rz(-1.3610164) q[2];
sx q[2];
rz(2.7854475) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1396461) q[1];
sx q[1];
rz(-1.0602385) q[1];
sx q[1];
rz(-0.21378784) q[1];
rz(-pi) q[2];
rz(1.5045937) q[3];
sx q[3];
rz(-2.8567064) q[3];
sx q[3];
rz(1.6411622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.88963738) q[2];
sx q[2];
rz(-0.95459443) q[2];
sx q[2];
rz(1.4944705) q[2];
rz(2.5877118) q[3];
sx q[3];
rz(-1.6051822) q[3];
sx q[3];
rz(1.4166098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81903356) q[0];
sx q[0];
rz(-0.77545866) q[0];
sx q[0];
rz(1.092528) q[0];
rz(-0.93387261) q[1];
sx q[1];
rz(-1.4051508) q[1];
sx q[1];
rz(0.78528231) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1500871) q[0];
sx q[0];
rz(-1.1508216) q[0];
sx q[0];
rz(0.49862592) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0243591) q[2];
sx q[2];
rz(-1.1866202) q[2];
sx q[2];
rz(0.4740997) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.45594564) q[1];
sx q[1];
rz(-0.37044493) q[1];
sx q[1];
rz(2.4042995) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5232072) q[3];
sx q[3];
rz(-1.9342285) q[3];
sx q[3];
rz(-1.5314764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2698722) q[2];
sx q[2];
rz(-1.4464658) q[2];
sx q[2];
rz(2.5093057) q[2];
rz(0.91442433) q[3];
sx q[3];
rz(-1.2605896) q[3];
sx q[3];
rz(-2.9818592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3226427) q[0];
sx q[0];
rz(-0.85700789) q[0];
sx q[0];
rz(-2.6981165) q[0];
rz(-0.30287287) q[1];
sx q[1];
rz(-2.5587406) q[1];
sx q[1];
rz(1.833896) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46638495) q[0];
sx q[0];
rz(-1.2699915) q[0];
sx q[0];
rz(1.3928901) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8820178) q[2];
sx q[2];
rz(-2.0577132) q[2];
sx q[2];
rz(-0.68070179) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.16567) q[1];
sx q[1];
rz(-2.5745086) q[1];
sx q[1];
rz(-2.0815064) q[1];
rz(-pi) q[2];
rz(-2.0409706) q[3];
sx q[3];
rz(-1.3047781) q[3];
sx q[3];
rz(3.0373552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.12094721) q[2];
sx q[2];
rz(-2.1006613) q[2];
sx q[2];
rz(-0.61332235) q[2];
rz(-0.81765085) q[3];
sx q[3];
rz(-0.48912564) q[3];
sx q[3];
rz(-0.31942719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15449512) q[0];
sx q[0];
rz(-1.7648062) q[0];
sx q[0];
rz(2.7823271) q[0];
rz(0.66419762) q[1];
sx q[1];
rz(-1.3281053) q[1];
sx q[1];
rz(0.42116234) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065124113) q[0];
sx q[0];
rz(-1.2026005) q[0];
sx q[0];
rz(0.23875956) q[0];
x q[1];
rz(3.1362651) q[2];
sx q[2];
rz(-1.6432646) q[2];
sx q[2];
rz(2.7203512) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2639271) q[1];
sx q[1];
rz(-1.3093728) q[1];
sx q[1];
rz(-1.1007453) q[1];
rz(-pi) q[2];
rz(-2.1499885) q[3];
sx q[3];
rz(-1.0458071) q[3];
sx q[3];
rz(0.92574173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26587129) q[2];
sx q[2];
rz(-3.044812) q[2];
sx q[2];
rz(-1.2145915) q[2];
rz(2.752979) q[3];
sx q[3];
rz(-0.92258421) q[3];
sx q[3];
rz(2.0428366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6471967) q[0];
sx q[0];
rz(-1.5955441) q[0];
sx q[0];
rz(1.2936976) q[0];
rz(2.9248059) q[1];
sx q[1];
rz(-0.68991359) q[1];
sx q[1];
rz(0.51295113) q[1];
rz(-2.4642261) q[2];
sx q[2];
rz(-1.2268885) q[2];
sx q[2];
rz(2.3033768) q[2];
rz(-2.67542) q[3];
sx q[3];
rz(-1.8099347) q[3];
sx q[3];
rz(3.1113887) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
