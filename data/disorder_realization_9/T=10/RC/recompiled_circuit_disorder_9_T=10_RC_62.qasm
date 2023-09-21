OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7115241) q[0];
sx q[0];
rz(-0.067458955) q[0];
sx q[0];
rz(0.67396069) q[0];
rz(2.824488) q[1];
sx q[1];
rz(-1.5082521) q[1];
sx q[1];
rz(-2.34692) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.062382467) q[0];
sx q[0];
rz(-1.198472) q[0];
sx q[0];
rz(1.1862434) q[0];
x q[1];
rz(1.7180213) q[2];
sx q[2];
rz(-1.0861673) q[2];
sx q[2];
rz(2.0369612) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4218581) q[1];
sx q[1];
rz(-1.0654447) q[1];
sx q[1];
rz(1.964633) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0178954) q[3];
sx q[3];
rz(-1.6959794) q[3];
sx q[3];
rz(-0.048429117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.47544605) q[2];
sx q[2];
rz(-2.011516) q[2];
sx q[2];
rz(-1.8165992) q[2];
rz(-2.825286) q[3];
sx q[3];
rz(-2.2241212) q[3];
sx q[3];
rz(0.90707183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9544202) q[0];
sx q[0];
rz(-0.81048727) q[0];
sx q[0];
rz(-0.95970884) q[0];
rz(-1.487544) q[1];
sx q[1];
rz(-0.68351775) q[1];
sx q[1];
rz(0.62746343) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42694416) q[0];
sx q[0];
rz(-1.5776002) q[0];
sx q[0];
rz(0.0077008458) q[0];
rz(-pi) q[1];
rz(0.10208315) q[2];
sx q[2];
rz(-0.21271579) q[2];
sx q[2];
rz(-1.4174457) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.56541601) q[1];
sx q[1];
rz(-1.7130573) q[1];
sx q[1];
rz(-1.8131282) q[1];
x q[2];
rz(-2.0262358) q[3];
sx q[3];
rz(-3.1133828) q[3];
sx q[3];
rz(-1.7445604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6015357) q[2];
sx q[2];
rz(-1.4011708) q[2];
sx q[2];
rz(0.97989782) q[2];
rz(1.4510138) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(-1.6987945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-0.56101218) q[0];
sx q[0];
rz(-3.09364) q[0];
sx q[0];
rz(-1.7618435) q[0];
rz(-0.17669949) q[1];
sx q[1];
rz(-1.9243762) q[1];
sx q[1];
rz(1.3476936) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4189258) q[0];
sx q[0];
rz(-1.330266) q[0];
sx q[0];
rz(0.47970432) q[0];
rz(-3.0323896) q[2];
sx q[2];
rz(-0.20143992) q[2];
sx q[2];
rz(-2.2815454) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8417931) q[1];
sx q[1];
rz(-1.6922957) q[1];
sx q[1];
rz(-1.6714765) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.090431902) q[3];
sx q[3];
rz(-1.7548429) q[3];
sx q[3];
rz(2.4726766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1220864) q[2];
sx q[2];
rz(-0.85947376) q[2];
sx q[2];
rz(-1.9879831) q[2];
rz(-2.7518318) q[3];
sx q[3];
rz(-0.4699769) q[3];
sx q[3];
rz(-0.78891689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6692114) q[0];
sx q[0];
rz(-2.5166469) q[0];
sx q[0];
rz(-0.56458449) q[0];
rz(3.0584884) q[1];
sx q[1];
rz(-1.9080947) q[1];
sx q[1];
rz(-0.21534236) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0523473) q[0];
sx q[0];
rz(-2.2627292) q[0];
sx q[0];
rz(-2.7034876) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16673659) q[2];
sx q[2];
rz(-0.98966375) q[2];
sx q[2];
rz(0.8596479) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.035186471) q[1];
sx q[1];
rz(-1.3760929) q[1];
sx q[1];
rz(0.10722843) q[1];
x q[2];
rz(1.1057042) q[3];
sx q[3];
rz(-1.0700738) q[3];
sx q[3];
rz(-2.6021258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5061491) q[2];
sx q[2];
rz(-2.052867) q[2];
sx q[2];
rz(-2.5396458) q[2];
rz(0.69058949) q[3];
sx q[3];
rz(-2.3875321) q[3];
sx q[3];
rz(-0.62600342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85201207) q[0];
sx q[0];
rz(-0.99442011) q[0];
sx q[0];
rz(2.1767298) q[0];
rz(-3.124974) q[1];
sx q[1];
rz(-0.42533541) q[1];
sx q[1];
rz(2.4749277) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0913038) q[0];
sx q[0];
rz(-1.6717523) q[0];
sx q[0];
rz(1.1703277) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7890116) q[2];
sx q[2];
rz(-2.0902299) q[2];
sx q[2];
rz(0.40371603) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.9469229) q[1];
sx q[1];
rz(-1.2877081) q[1];
sx q[1];
rz(-2.4592295) q[1];
x q[2];
rz(2.9186451) q[3];
sx q[3];
rz(-1.8604391) q[3];
sx q[3];
rz(1.2900316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2255286) q[2];
sx q[2];
rz(-1.7614044) q[2];
sx q[2];
rz(1.2529681) q[2];
rz(1.3145674) q[3];
sx q[3];
rz(-1.9606749) q[3];
sx q[3];
rz(-2.0107646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1154293) q[0];
sx q[0];
rz(-2.207343) q[0];
sx q[0];
rz(-1.3622267) q[0];
rz(-1.7419787) q[1];
sx q[1];
rz(-1.5618961) q[1];
sx q[1];
rz(-1.5302352) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0936733) q[0];
sx q[0];
rz(-1.9903272) q[0];
sx q[0];
rz(1.6506667) q[0];
rz(-pi) q[1];
rz(1.2457232) q[2];
sx q[2];
rz(-1.1522066) q[2];
sx q[2];
rz(-0.77667728) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6758319) q[1];
sx q[1];
rz(-2.6869046) q[1];
sx q[1];
rz(2.0783706) q[1];
rz(-pi) q[2];
rz(-0.050614428) q[3];
sx q[3];
rz(-1.317987) q[3];
sx q[3];
rz(1.3150172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.35649148) q[2];
sx q[2];
rz(-1.0078398) q[2];
sx q[2];
rz(0.029416857) q[2];
rz(1.6507089) q[3];
sx q[3];
rz(-2.1865032) q[3];
sx q[3];
rz(1.9194549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71762639) q[0];
sx q[0];
rz(-0.76537776) q[0];
sx q[0];
rz(-0.18280612) q[0];
rz(-0.43552661) q[1];
sx q[1];
rz(-0.85414779) q[1];
sx q[1];
rz(2.6307154) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8144007) q[0];
sx q[0];
rz(-1.4609903) q[0];
sx q[0];
rz(0.26392428) q[0];
rz(-2.3350231) q[2];
sx q[2];
rz(-0.91420805) q[2];
sx q[2];
rz(-0.29733959) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6888914) q[1];
sx q[1];
rz(-2.211262) q[1];
sx q[1];
rz(-1.0673317) q[1];
rz(-pi) q[2];
rz(2.0914518) q[3];
sx q[3];
rz(-2.5556892) q[3];
sx q[3];
rz(1.0019765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9110979) q[2];
sx q[2];
rz(-0.7926597) q[2];
sx q[2];
rz(-0.62057173) q[2];
rz(-0.42256045) q[3];
sx q[3];
rz(-1.3603323) q[3];
sx q[3];
rz(1.871199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34185103) q[0];
sx q[0];
rz(-0.36062476) q[0];
sx q[0];
rz(-2.0525232) q[0];
rz(-1.0808806) q[1];
sx q[1];
rz(-1.2842478) q[1];
sx q[1];
rz(0.63527766) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7437744) q[0];
sx q[0];
rz(-1.213307) q[0];
sx q[0];
rz(2.9767569) q[0];
x q[1];
rz(-0.17417553) q[2];
sx q[2];
rz(-2.0617699) q[2];
sx q[2];
rz(-0.27506405) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2891846) q[1];
sx q[1];
rz(-1.8361366) q[1];
sx q[1];
rz(-0.62383382) q[1];
rz(0.99526309) q[3];
sx q[3];
rz(-1.4008153) q[3];
sx q[3];
rz(0.02804027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0044331) q[2];
sx q[2];
rz(-1.6798423) q[2];
sx q[2];
rz(-3.1414202) q[2];
rz(-2.4553283) q[3];
sx q[3];
rz(-0.72824794) q[3];
sx q[3];
rz(0.85787684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8453318) q[0];
sx q[0];
rz(-0.91264549) q[0];
sx q[0];
rz(2.4849179) q[0];
rz(2.7776921) q[1];
sx q[1];
rz(-0.44635043) q[1];
sx q[1];
rz(0.90973967) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0168403) q[0];
sx q[0];
rz(-2.7144055) q[0];
sx q[0];
rz(-2.4401526) q[0];
rz(-pi) q[1];
rz(1.3526731) q[2];
sx q[2];
rz(-1.2303196) q[2];
sx q[2];
rz(3.043963) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4917131) q[1];
sx q[1];
rz(-2.8976106) q[1];
sx q[1];
rz(2.6073543) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2551366) q[3];
sx q[3];
rz(-1.8161895) q[3];
sx q[3];
rz(2.0592225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3108814) q[2];
sx q[2];
rz(-2.723366) q[2];
sx q[2];
rz(0.48661423) q[2];
rz(3.0330372) q[3];
sx q[3];
rz(-2.4311549) q[3];
sx q[3];
rz(-1.2238812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20031032) q[0];
sx q[0];
rz(-2.0599984) q[0];
sx q[0];
rz(0.51668984) q[0];
rz(1.5962881) q[1];
sx q[1];
rz(-2.089112) q[1];
sx q[1];
rz(-2.2470078) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6870118) q[0];
sx q[0];
rz(-0.63475383) q[0];
sx q[0];
rz(1.6420341) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1799699) q[2];
sx q[2];
rz(-0.50190364) q[2];
sx q[2];
rz(-0.48197907) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.91720944) q[1];
sx q[1];
rz(-1.50748) q[1];
sx q[1];
rz(-1.9320095) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9856332) q[3];
sx q[3];
rz(-2.2150063) q[3];
sx q[3];
rz(-0.19176602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9937925) q[2];
sx q[2];
rz(-2.1788969) q[2];
sx q[2];
rz(2.1195892) q[2];
rz(1.3379124) q[3];
sx q[3];
rz(-1.6513848) q[3];
sx q[3];
rz(0.62906229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5932896) q[0];
sx q[0];
rz(-1.6175445) q[0];
sx q[0];
rz(2.2486726) q[0];
rz(-1.0409566) q[1];
sx q[1];
rz(-3.0491842) q[1];
sx q[1];
rz(-1.4728117) q[1];
rz(-2.5953318) q[2];
sx q[2];
rz(-0.98304521) q[2];
sx q[2];
rz(1.8084768) q[2];
rz(1.8716807) q[3];
sx q[3];
rz(-1.6172998) q[3];
sx q[3];
rz(-0.87340364) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
