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
rz(0.25843698) q[0];
sx q[0];
rz(-0.28764495) q[0];
sx q[0];
rz(-1.1076466) q[0];
rz(-1.8183174) q[1];
sx q[1];
rz(-0.32185093) q[1];
sx q[1];
rz(-2.5133361) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4391651) q[0];
sx q[0];
rz(-1.747923) q[0];
sx q[0];
rz(-0.11985531) q[0];
rz(-pi) q[1];
x q[1];
rz(0.722856) q[2];
sx q[2];
rz(-1.8713163) q[2];
sx q[2];
rz(2.3006257) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3648277) q[1];
sx q[1];
rz(-1.3421632) q[1];
sx q[1];
rz(1.6079812) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2443549) q[3];
sx q[3];
rz(-2.4938886) q[3];
sx q[3];
rz(3.0239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2884752) q[2];
sx q[2];
rz(-1.4318117) q[2];
sx q[2];
rz(-2.8004069) q[2];
rz(2.2513385) q[3];
sx q[3];
rz(-1.2672) q[3];
sx q[3];
rz(2.6064579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61966908) q[0];
sx q[0];
rz(-2.4578019) q[0];
sx q[0];
rz(-2.4271915) q[0];
rz(-1.5996541) q[1];
sx q[1];
rz(-2.6456412) q[1];
sx q[1];
rz(-0.094706789) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45825645) q[0];
sx q[0];
rz(-0.33412513) q[0];
sx q[0];
rz(-1.279356) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96927283) q[2];
sx q[2];
rz(-1.8516632) q[2];
sx q[2];
rz(-1.2796677) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.73037877) q[1];
sx q[1];
rz(-1.8068562) q[1];
sx q[1];
rz(-2.7884931) q[1];
x q[2];
rz(0.27908947) q[3];
sx q[3];
rz(-2.6153339) q[3];
sx q[3];
rz(2.581439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9793205) q[2];
sx q[2];
rz(-1.2459735) q[2];
sx q[2];
rz(2.6675513) q[2];
rz(2.3794409) q[3];
sx q[3];
rz(-2.5049329) q[3];
sx q[3];
rz(-0.4711841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8629465) q[0];
sx q[0];
rz(-0.11676783) q[0];
sx q[0];
rz(0.85292029) q[0];
rz(-2.9184753) q[1];
sx q[1];
rz(-2.6631963) q[1];
sx q[1];
rz(-0.51582897) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4321126) q[0];
sx q[0];
rz(-1.6211544) q[0];
sx q[0];
rz(-1.5169282) q[0];
rz(0.93348301) q[2];
sx q[2];
rz(-1.3636944) q[2];
sx q[2];
rz(2.0014971) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.22166907) q[1];
sx q[1];
rz(-0.95429568) q[1];
sx q[1];
rz(2.4024582) q[1];
rz(-pi) q[2];
rz(-2.2688341) q[3];
sx q[3];
rz(-1.8379194) q[3];
sx q[3];
rz(-1.1760548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8082661) q[2];
sx q[2];
rz(-0.9767248) q[2];
sx q[2];
rz(3.0466363) q[2];
rz(2.700108) q[3];
sx q[3];
rz(-1.35651) q[3];
sx q[3];
rz(1.1536185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74986356) q[0];
sx q[0];
rz(-2.5666105) q[0];
sx q[0];
rz(-2.0854501) q[0];
rz(-0.9437584) q[1];
sx q[1];
rz(-0.25139233) q[1];
sx q[1];
rz(1.4518849) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7923342) q[0];
sx q[0];
rz(-0.99674388) q[0];
sx q[0];
rz(-0.7416772) q[0];
rz(-pi) q[1];
rz(-0.69047549) q[2];
sx q[2];
rz(-2.568141) q[2];
sx q[2];
rz(0.8586463) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.02921) q[1];
sx q[1];
rz(-1.60648) q[1];
sx q[1];
rz(0.24072631) q[1];
x q[2];
rz(3.0309276) q[3];
sx q[3];
rz(-1.5783211) q[3];
sx q[3];
rz(1.8845817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7283432) q[2];
sx q[2];
rz(-2.3409833) q[2];
sx q[2];
rz(2.441414) q[2];
rz(0.87107301) q[3];
sx q[3];
rz(-2.0756523) q[3];
sx q[3];
rz(1.7843436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0102608) q[0];
sx q[0];
rz(-2.6483674) q[0];
sx q[0];
rz(2.6569271) q[0];
rz(2.2635745) q[1];
sx q[1];
rz(-1.1596102) q[1];
sx q[1];
rz(-2.7427618) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3225375) q[0];
sx q[0];
rz(-2.5761671) q[0];
sx q[0];
rz(-0.25382332) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9100175) q[2];
sx q[2];
rz(-0.89618635) q[2];
sx q[2];
rz(2.2774027) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5113885) q[1];
sx q[1];
rz(-0.92899073) q[1];
sx q[1];
rz(2.6766269) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2696617) q[3];
sx q[3];
rz(-2.0060349) q[3];
sx q[3];
rz(-0.98952451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.47348076) q[2];
sx q[2];
rz(-1.6010189) q[2];
sx q[2];
rz(0.87781805) q[2];
rz(2.7569438) q[3];
sx q[3];
rz(-2.2967702) q[3];
sx q[3];
rz(1.1788684) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28410742) q[0];
sx q[0];
rz(-0.50898886) q[0];
sx q[0];
rz(-2.8711163) q[0];
rz(-0.30666223) q[1];
sx q[1];
rz(-0.51391927) q[1];
sx q[1];
rz(-2.564863) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33615112) q[0];
sx q[0];
rz(-0.2052983) q[0];
sx q[0];
rz(-2.0144256) q[0];
rz(-pi) q[1];
x q[1];
rz(2.001695) q[2];
sx q[2];
rz(-2.4043519) q[2];
sx q[2];
rz(-0.37079045) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.19837436) q[1];
sx q[1];
rz(-0.70462185) q[1];
sx q[1];
rz(0.93438973) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33476013) q[3];
sx q[3];
rz(-1.5575512) q[3];
sx q[3];
rz(0.56559004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7755255) q[2];
sx q[2];
rz(-1.0696573) q[2];
sx q[2];
rz(2.1712187) q[2];
rz(2.7905285) q[3];
sx q[3];
rz(-0.23215663) q[3];
sx q[3];
rz(-2.8620913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2358667) q[0];
sx q[0];
rz(-0.38385639) q[0];
sx q[0];
rz(-0.50091499) q[0];
rz(-2.4210335) q[1];
sx q[1];
rz(-0.84545285) q[1];
sx q[1];
rz(-0.93773425) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7739627) q[0];
sx q[0];
rz(-1.367462) q[0];
sx q[0];
rz(-0.042026566) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9483614) q[2];
sx q[2];
rz(-2.015471) q[2];
sx q[2];
rz(-1.3832165) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5590458) q[1];
sx q[1];
rz(-1.2101047) q[1];
sx q[1];
rz(-2.2466895) q[1];
rz(3.1295007) q[3];
sx q[3];
rz(-2.4547057) q[3];
sx q[3];
rz(-1.6364607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0096036) q[2];
sx q[2];
rz(-0.14543532) q[2];
sx q[2];
rz(2.2930938) q[2];
rz(2.2961473) q[3];
sx q[3];
rz(-0.65631056) q[3];
sx q[3];
rz(1.2257303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7312412) q[0];
sx q[0];
rz(-0.9587962) q[0];
sx q[0];
rz(2.248726) q[0];
rz(-3.0800173) q[1];
sx q[1];
rz(-2.4696746) q[1];
sx q[1];
rz(-0.16199131) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0688871) q[0];
sx q[0];
rz(-1.6729682) q[0];
sx q[0];
rz(-3.1260558) q[0];
rz(-pi) q[1];
rz(-1.685437) q[2];
sx q[2];
rz(-0.87213665) q[2];
sx q[2];
rz(1.4574775) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3969035) q[1];
sx q[1];
rz(-1.2619881) q[1];
sx q[1];
rz(-1.7148605) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4668767) q[3];
sx q[3];
rz(-1.4774408) q[3];
sx q[3];
rz(-2.4568077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.050345) q[2];
sx q[2];
rz(-0.42751905) q[2];
sx q[2];
rz(0.72489911) q[2];
rz(-0.27724087) q[3];
sx q[3];
rz(-0.96479779) q[3];
sx q[3];
rz(-3.0550756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25170931) q[0];
sx q[0];
rz(-0.6530264) q[0];
sx q[0];
rz(3.1157893) q[0];
rz(1.7426527) q[1];
sx q[1];
rz(-1.502123) q[1];
sx q[1];
rz(0.10765156) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9190034) q[0];
sx q[0];
rz(-2.1371142) q[0];
sx q[0];
rz(-2.983271) q[0];
rz(0.74371643) q[2];
sx q[2];
rz(-2.5355843) q[2];
sx q[2];
rz(-1.5845989) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4835158) q[1];
sx q[1];
rz(-0.47464322) q[1];
sx q[1];
rz(-0.078981699) q[1];
x q[2];
rz(-0.017849542) q[3];
sx q[3];
rz(-1.4571725) q[3];
sx q[3];
rz(-2.2872432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.671635) q[2];
sx q[2];
rz(-1.8822972) q[2];
sx q[2];
rz(-1.0939481) q[2];
rz(0.57540244) q[3];
sx q[3];
rz(-2.3016774) q[3];
sx q[3];
rz(1.1168787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53330082) q[0];
sx q[0];
rz(-0.58654439) q[0];
sx q[0];
rz(-2.4285512) q[0];
rz(0.22480045) q[1];
sx q[1];
rz(-2.5495922) q[1];
sx q[1];
rz(-1.0932659) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9295392) q[0];
sx q[0];
rz(-2.0418379) q[0];
sx q[0];
rz(-2.3624198) q[0];
rz(2.9399038) q[2];
sx q[2];
rz(-1.301577) q[2];
sx q[2];
rz(-1.2303011) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.374129) q[1];
sx q[1];
rz(-1.1157562) q[1];
sx q[1];
rz(-1.6947075) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8328463) q[3];
sx q[3];
rz(-2.7152677) q[3];
sx q[3];
rz(-2.654512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1405868) q[2];
sx q[2];
rz(-0.69585496) q[2];
sx q[2];
rz(-0.33983964) q[2];
rz(0.042424399) q[3];
sx q[3];
rz(-1.2845311) q[3];
sx q[3];
rz(-2.5137918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9545659) q[0];
sx q[0];
rz(-1.621959) q[0];
sx q[0];
rz(1.8931615) q[0];
rz(-2.3595702) q[1];
sx q[1];
rz(-2.2047058) q[1];
sx q[1];
rz(2.4155736) q[1];
rz(0.37219477) q[2];
sx q[2];
rz(-1.8831461) q[2];
sx q[2];
rz(2.1103721) q[2];
rz(-1.9966765) q[3];
sx q[3];
rz(-1.799634) q[3];
sx q[3];
rz(-1.4992731) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
