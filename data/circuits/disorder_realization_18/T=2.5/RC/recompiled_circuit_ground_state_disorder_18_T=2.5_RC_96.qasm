OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6567) q[0];
sx q[0];
rz(-0.088887366) q[0];
sx q[0];
rz(2.3645997) q[0];
rz(0.69500336) q[1];
sx q[1];
rz(-3.0106795) q[1];
sx q[1];
rz(-2.793269) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3366813) q[0];
sx q[0];
rz(-0.8136734) q[0];
sx q[0];
rz(1.8437852) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.25599249) q[2];
sx q[2];
rz(-2.8446994) q[2];
sx q[2];
rz(2.1741637) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4433044) q[1];
sx q[1];
rz(-2.0870247) q[1];
sx q[1];
rz(-0.97133941) q[1];
rz(-pi) q[2];
rz(-1.6197085) q[3];
sx q[3];
rz(-1.1756611) q[3];
sx q[3];
rz(1.6312158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6473039) q[2];
sx q[2];
rz(-0.56607294) q[2];
sx q[2];
rz(-1.0757793) q[2];
rz(3.0268269) q[3];
sx q[3];
rz(-1.4817295) q[3];
sx q[3];
rz(2.4858294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7455416) q[0];
sx q[0];
rz(-2.2069554) q[0];
sx q[0];
rz(-2.4387687) q[0];
rz(-1.5463691) q[1];
sx q[1];
rz(-0.74235761) q[1];
sx q[1];
rz(-2.5667618) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2947161) q[0];
sx q[0];
rz(-2.623577) q[0];
sx q[0];
rz(-2.4344492) q[0];
x q[1];
rz(2.4509239) q[2];
sx q[2];
rz(-2.4021742) q[2];
sx q[2];
rz(-1.2659466) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.489255) q[1];
sx q[1];
rz(-1.7018507) q[1];
sx q[1];
rz(-0.62789692) q[1];
x q[2];
rz(-2.094926) q[3];
sx q[3];
rz(-1.0051749) q[3];
sx q[3];
rz(-0.079295955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.90291643) q[2];
sx q[2];
rz(-1.6409589) q[2];
sx q[2];
rz(-2.404876) q[2];
rz(-2.369407) q[3];
sx q[3];
rz(-0.18852791) q[3];
sx q[3];
rz(-2.7948715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(2.9521028) q[0];
sx q[0];
rz(-1.075407) q[0];
sx q[0];
rz(1.1601675) q[0];
rz(0.53933764) q[1];
sx q[1];
rz(-0.62349206) q[1];
sx q[1];
rz(-3.0711874) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6269913) q[0];
sx q[0];
rz(-0.87341269) q[0];
sx q[0];
rz(-2.227289) q[0];
rz(-pi) q[1];
rz(2.3023476) q[2];
sx q[2];
rz(-2.9021724) q[2];
sx q[2];
rz(-0.024062238) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.119795) q[1];
sx q[1];
rz(-2.7547495) q[1];
sx q[1];
rz(-2.9683044) q[1];
x q[2];
rz(-1.0892404) q[3];
sx q[3];
rz(-2.1641971) q[3];
sx q[3];
rz(-2.7508222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.77014852) q[2];
sx q[2];
rz(-1.4263209) q[2];
sx q[2];
rz(-2.9798129) q[2];
rz(-0.77156228) q[3];
sx q[3];
rz(-3.0050889) q[3];
sx q[3];
rz(2.0989044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55018392) q[0];
sx q[0];
rz(-2.3593481) q[0];
sx q[0];
rz(-0.25006205) q[0];
rz(1.0391883) q[1];
sx q[1];
rz(-2.3999374) q[1];
sx q[1];
rz(0.80113062) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.259819) q[0];
sx q[0];
rz(-1.7823813) q[0];
sx q[0];
rz(-1.5552646) q[0];
rz(-pi) q[1];
rz(-0.25867985) q[2];
sx q[2];
rz(-1.341267) q[2];
sx q[2];
rz(1.2327641) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.86601195) q[1];
sx q[1];
rz(-0.51803127) q[1];
sx q[1];
rz(-0.14054246) q[1];
rz(-0.43901033) q[3];
sx q[3];
rz(-1.9042104) q[3];
sx q[3];
rz(-0.96102137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.73418009) q[2];
sx q[2];
rz(-1.696618) q[2];
sx q[2];
rz(2.5511191) q[2];
rz(-0.34481314) q[3];
sx q[3];
rz(-2.7761288) q[3];
sx q[3];
rz(-1.6706246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.19584) q[0];
sx q[0];
rz(-1.4386289) q[0];
sx q[0];
rz(1.6666743) q[0];
rz(-1.6189812) q[1];
sx q[1];
rz(-1.5432065) q[1];
sx q[1];
rz(-2.7326857) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4034861) q[0];
sx q[0];
rz(-0.047560725) q[0];
sx q[0];
rz(-0.17769583) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7402503) q[2];
sx q[2];
rz(-2.231488) q[2];
sx q[2];
rz(-1.5398538) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6521367) q[1];
sx q[1];
rz(-2.5083275) q[1];
sx q[1];
rz(-0.79852028) q[1];
rz(-pi) q[2];
rz(-1.4765396) q[3];
sx q[3];
rz(-0.89897147) q[3];
sx q[3];
rz(-2.6420324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.54395479) q[2];
sx q[2];
rz(-3.0089617) q[2];
sx q[2];
rz(-0.65072101) q[2];
rz(1.4607653) q[3];
sx q[3];
rz(-1.4131578) q[3];
sx q[3];
rz(2.7248603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2862709) q[0];
sx q[0];
rz(-2.0725508) q[0];
sx q[0];
rz(-1.8977813) q[0];
rz(0.52737251) q[1];
sx q[1];
rz(-1.6970044) q[1];
sx q[1];
rz(1.3810371) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6207989) q[0];
sx q[0];
rz(-0.24220151) q[0];
sx q[0];
rz(0.050817055) q[0];
rz(1.6445134) q[2];
sx q[2];
rz(-1.4716765) q[2];
sx q[2];
rz(0.69621554) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.11070793) q[1];
sx q[1];
rz(-1.0652911) q[1];
sx q[1];
rz(0.55191374) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1072371) q[3];
sx q[3];
rz(-1.966668) q[3];
sx q[3];
rz(-0.27908868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0397628) q[2];
sx q[2];
rz(-1.3527752) q[2];
sx q[2];
rz(-1.1438032) q[2];
rz(2.3666798) q[3];
sx q[3];
rz(-1.9655922) q[3];
sx q[3];
rz(2.9162143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5932805) q[0];
sx q[0];
rz(-0.53638023) q[0];
sx q[0];
rz(-0.29658741) q[0];
rz(-2.2071154) q[1];
sx q[1];
rz(-2.2521033) q[1];
sx q[1];
rz(2.8737822) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4992384) q[0];
sx q[0];
rz(-1.7846414) q[0];
sx q[0];
rz(1.3539899) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8765062) q[2];
sx q[2];
rz(-2.5137182) q[2];
sx q[2];
rz(-1.6123079) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0632656) q[1];
sx q[1];
rz(-0.73042471) q[1];
sx q[1];
rz(0.46565454) q[1];
rz(0.15800627) q[3];
sx q[3];
rz(-1.7078764) q[3];
sx q[3];
rz(-0.44120818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3725793) q[2];
sx q[2];
rz(-2.09435) q[2];
sx q[2];
rz(-3.0229342) q[2];
rz(-0.81698051) q[3];
sx q[3];
rz(-1.1840362) q[3];
sx q[3];
rz(0.55678862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34767041) q[0];
sx q[0];
rz(-0.477328) q[0];
sx q[0];
rz(-1.073786) q[0];
rz(-2.1444164) q[1];
sx q[1];
rz(-1.8938096) q[1];
sx q[1];
rz(1.006385) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2074315) q[0];
sx q[0];
rz(-1.6994358) q[0];
sx q[0];
rz(-3.1324803) q[0];
x q[1];
rz(-0.7147737) q[2];
sx q[2];
rz(-2.6539485) q[2];
sx q[2];
rz(-0.89370103) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.83676618) q[1];
sx q[1];
rz(-1.8988892) q[1];
sx q[1];
rz(1.3681218) q[1];
rz(0.96285497) q[3];
sx q[3];
rz(-0.7391151) q[3];
sx q[3];
rz(-1.407042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5239079) q[2];
sx q[2];
rz(-2.0670117) q[2];
sx q[2];
rz(-2.1802444) q[2];
rz(-2.5210467) q[3];
sx q[3];
rz(-2.4397662) q[3];
sx q[3];
rz(0.8121686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0470225) q[0];
sx q[0];
rz(-0.71806878) q[0];
sx q[0];
rz(0.72599232) q[0];
rz(-0.30966169) q[1];
sx q[1];
rz(-1.0331215) q[1];
sx q[1];
rz(2.9594701) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5512729) q[0];
sx q[0];
rz(-2.3819551) q[0];
sx q[0];
rz(2.9653427) q[0];
x q[1];
rz(-2.0810037) q[2];
sx q[2];
rz(-1.1119208) q[2];
sx q[2];
rz(-1.4993678) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9023064) q[1];
sx q[1];
rz(-0.42879802) q[1];
sx q[1];
rz(2.8132755) q[1];
rz(-0.52892367) q[3];
sx q[3];
rz(-0.70883239) q[3];
sx q[3];
rz(0.20656997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9228268) q[2];
sx q[2];
rz(-1.7106235) q[2];
sx q[2];
rz(1.1440841) q[2];
rz(0.73793441) q[3];
sx q[3];
rz(-2.4149371) q[3];
sx q[3];
rz(-1.1206333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4436307) q[0];
sx q[0];
rz(-1.8101036) q[0];
sx q[0];
rz(-0.6829845) q[0];
rz(1.5890652) q[1];
sx q[1];
rz(-2.2582159) q[1];
sx q[1];
rz(2.24486) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6370332) q[0];
sx q[0];
rz(-2.6650223) q[0];
sx q[0];
rz(-0.43463114) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.801034) q[2];
sx q[2];
rz(-0.55998324) q[2];
sx q[2];
rz(-0.79027938) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7234133) q[1];
sx q[1];
rz(-1.6589983) q[1];
sx q[1];
rz(-3.0107016) q[1];
rz(-pi) q[2];
rz(2.0192782) q[3];
sx q[3];
rz(-2.0297673) q[3];
sx q[3];
rz(-2.36623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5985742) q[2];
sx q[2];
rz(-1.2335346) q[2];
sx q[2];
rz(-0.42178556) q[2];
rz(-0.56454033) q[3];
sx q[3];
rz(-0.83344236) q[3];
sx q[3];
rz(2.25901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-0.17557872) q[0];
sx q[0];
rz(-1.6015552) q[0];
sx q[0];
rz(1.8764499) q[0];
rz(-1.235678) q[1];
sx q[1];
rz(-1.0503678) q[1];
sx q[1];
rz(-2.8138524) q[1];
rz(2.7719978) q[2];
sx q[2];
rz(-1.3736178) q[2];
sx q[2];
rz(-1.553229) q[2];
rz(2.3269736) q[3];
sx q[3];
rz(-0.64523166) q[3];
sx q[3];
rz(1.393691) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
