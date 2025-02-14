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
rz(2.0339461) q[0];
rz(1.3232752) q[1];
sx q[1];
rz(-2.8197417) q[1];
sx q[1];
rz(-0.62825656) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88958607) q[0];
sx q[0];
rz(-1.4528251) q[0];
sx q[0];
rz(-1.3924166) q[0];
x q[1];
rz(0.43810644) q[2];
sx q[2];
rz(-0.77229653) q[2];
sx q[2];
rz(-0.4060678) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7767649) q[1];
sx q[1];
rz(-1.3421632) q[1];
sx q[1];
rz(-1.5336114) q[1];
x q[2];
rz(-2.2443549) q[3];
sx q[3];
rz(-2.4938886) q[3];
sx q[3];
rz(-0.11769262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2884752) q[2];
sx q[2];
rz(-1.709781) q[2];
sx q[2];
rz(-2.8004069) q[2];
rz(2.2513385) q[3];
sx q[3];
rz(-1.8743926) q[3];
sx q[3];
rz(0.53513479) q[3];
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
rz(-pi) q[0];
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
rz(-0.61966908) q[0];
sx q[0];
rz(-0.6837908) q[0];
sx q[0];
rz(2.4271915) q[0];
rz(1.5996541) q[1];
sx q[1];
rz(-2.6456412) q[1];
sx q[1];
rz(-3.0468859) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7529012) q[0];
sx q[0];
rz(-1.6651648) q[0];
sx q[0];
rz(-1.8917985) q[0];
rz(0.33659597) q[2];
sx q[2];
rz(-2.1456652) q[2];
sx q[2];
rz(-0.47913715) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.73037877) q[1];
sx q[1];
rz(-1.8068562) q[1];
sx q[1];
rz(2.7884931) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4121145) q[3];
sx q[3];
rz(-2.074721) q[3];
sx q[3];
rz(-0.24013396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16227214) q[2];
sx q[2];
rz(-1.8956192) q[2];
sx q[2];
rz(2.6675513) q[2];
rz(-0.76215172) q[3];
sx q[3];
rz(-2.5049329) q[3];
sx q[3];
rz(-0.4711841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8629465) q[0];
sx q[0];
rz(-0.11676783) q[0];
sx q[0];
rz(-2.2886724) q[0];
rz(-2.9184753) q[1];
sx q[1];
rz(-2.6631963) q[1];
sx q[1];
rz(-0.51582897) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2775622) q[0];
sx q[0];
rz(-1.5169965) q[0];
sx q[0];
rz(0.050431099) q[0];
rz(1.91024) q[2];
sx q[2];
rz(-2.4759549) q[2];
sx q[2];
rz(-2.9817944) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9199236) q[1];
sx q[1];
rz(-0.95429568) q[1];
sx q[1];
rz(-0.73913445) q[1];
x q[2];
rz(-0.34308691) q[3];
sx q[3];
rz(-2.2394387) q[3];
sx q[3];
rz(2.9647788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3333266) q[2];
sx q[2];
rz(-2.1648679) q[2];
sx q[2];
rz(0.094956368) q[2];
rz(0.44148463) q[3];
sx q[3];
rz(-1.7850826) q[3];
sx q[3];
rz(1.1536185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74986356) q[0];
sx q[0];
rz(-0.57498217) q[0];
sx q[0];
rz(1.0561426) q[0];
rz(0.9437584) q[1];
sx q[1];
rz(-0.25139233) q[1];
sx q[1];
rz(-1.4518849) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3817134) q[0];
sx q[0];
rz(-0.96764123) q[0];
sx q[0];
rz(-2.2908014) q[0];
rz(1.9610434) q[2];
sx q[2];
rz(-2.0023291) q[2];
sx q[2];
rz(0.081588946) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8274967) q[1];
sx q[1];
rz(-2.8982867) q[1];
sx q[1];
rz(-0.14863405) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5632252) q[3];
sx q[3];
rz(-1.6814582) q[3];
sx q[3];
rz(-0.31294926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4132495) q[2];
sx q[2];
rz(-2.3409833) q[2];
sx q[2];
rz(-2.441414) q[2];
rz(-0.87107301) q[3];
sx q[3];
rz(-2.0756523) q[3];
sx q[3];
rz(1.357249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1313318) q[0];
sx q[0];
rz(-0.49322525) q[0];
sx q[0];
rz(0.48466551) q[0];
rz(2.2635745) q[1];
sx q[1];
rz(-1.9819825) q[1];
sx q[1];
rz(-0.39883089) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53610401) q[0];
sx q[0];
rz(-1.4358504) q[0];
sx q[0];
rz(0.55079726) q[0];
rz(1.9100175) q[2];
sx q[2];
rz(-2.2454063) q[2];
sx q[2];
rz(-0.86418992) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64885253) q[1];
sx q[1];
rz(-1.938174) q[1];
sx q[1];
rz(0.87441625) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5957803) q[3];
sx q[3];
rz(-0.94792507) q[3];
sx q[3];
rz(0.92178492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.47348076) q[2];
sx q[2];
rz(-1.5405737) q[2];
sx q[2];
rz(-0.87781805) q[2];
rz(-0.38464883) q[3];
sx q[3];
rz(-0.84482241) q[3];
sx q[3];
rz(-1.1788684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28410742) q[0];
sx q[0];
rz(-0.50898886) q[0];
sx q[0];
rz(2.8711163) q[0];
rz(2.8349304) q[1];
sx q[1];
rz(-2.6276734) q[1];
sx q[1];
rz(-0.57672966) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33615112) q[0];
sx q[0];
rz(-2.9362943) q[0];
sx q[0];
rz(2.0144256) q[0];
rz(-pi) q[1];
rz(-2.001695) q[2];
sx q[2];
rz(-0.73724079) q[2];
sx q[2];
rz(-0.37079045) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.19837436) q[1];
sx q[1];
rz(-0.70462185) q[1];
sx q[1];
rz(0.93438973) q[1];
rz(-pi) q[2];
rz(-2.8068325) q[3];
sx q[3];
rz(-1.5840414) q[3];
sx q[3];
rz(2.5760026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3660672) q[2];
sx q[2];
rz(-1.0696573) q[2];
sx q[2];
rz(2.1712187) q[2];
rz(0.35106418) q[3];
sx q[3];
rz(-0.23215663) q[3];
sx q[3];
rz(2.8620913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2358667) q[0];
sx q[0];
rz(-2.7577363) q[0];
sx q[0];
rz(-2.6406777) q[0];
rz(-2.4210335) q[1];
sx q[1];
rz(-2.2961398) q[1];
sx q[1];
rz(-2.2038584) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7739627) q[0];
sx q[0];
rz(-1.7741307) q[0];
sx q[0];
rz(0.042026566) q[0];
x q[1];
rz(-2.6679084) q[2];
sx q[2];
rz(-1.9100744) q[2];
sx q[2];
rz(2.7850399) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.42635065) q[1];
sx q[1];
rz(-2.3890426) q[1];
sx q[1];
rz(2.1133711) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68685112) q[3];
sx q[3];
rz(-1.5631286) q[3];
sx q[3];
rz(0.075014278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0096036) q[2];
sx q[2];
rz(-2.9961573) q[2];
sx q[2];
rz(2.2930938) q[2];
rz(0.84544539) q[3];
sx q[3];
rz(-2.4852821) q[3];
sx q[3];
rz(-1.9158624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41035143) q[0];
sx q[0];
rz(-0.9587962) q[0];
sx q[0];
rz(-2.248726) q[0];
rz(-3.0800173) q[1];
sx q[1];
rz(-2.4696746) q[1];
sx q[1];
rz(-0.16199131) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2238857) q[0];
sx q[0];
rz(-3.0382503) q[0];
sx q[0];
rz(-1.4204106) q[0];
rz(-pi) q[1];
rz(-3.0062469) q[2];
sx q[2];
rz(-0.70643808) q[2];
sx q[2];
rz(1.634623) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8422441) q[1];
sx q[1];
rz(-0.33978251) q[1];
sx q[1];
rz(-2.7187125) q[1];
rz(-pi) q[2];
rz(-2.3050383) q[3];
sx q[3];
rz(-0.13958195) q[3];
sx q[3];
rz(0.15662801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.091247678) q[2];
sx q[2];
rz(-2.7140736) q[2];
sx q[2];
rz(-0.72489911) q[2];
rz(-2.8643518) q[3];
sx q[3];
rz(-2.1767949) q[3];
sx q[3];
rz(0.086517081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25170931) q[0];
sx q[0];
rz(-0.6530264) q[0];
sx q[0];
rz(0.025803331) q[0];
rz(1.39894) q[1];
sx q[1];
rz(-1.502123) q[1];
sx q[1];
rz(3.0339411) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9190034) q[0];
sx q[0];
rz(-2.1371142) q[0];
sx q[0];
rz(0.1583216) q[0];
x q[1];
rz(-1.1321105) q[2];
sx q[2];
rz(-2.0033547) q[2];
sx q[2];
rz(0.74299845) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.15756179) q[1];
sx q[1];
rz(-1.6068629) q[1];
sx q[1];
rz(-0.47337516) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4571545) q[3];
sx q[3];
rz(-1.5530619) q[3];
sx q[3];
rz(2.4271698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.671635) q[2];
sx q[2];
rz(-1.8822972) q[2];
sx q[2];
rz(-2.0476445) q[2];
rz(-2.5661902) q[3];
sx q[3];
rz(-0.83991528) q[3];
sx q[3];
rz(2.0247139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53330082) q[0];
sx q[0];
rz(-0.58654439) q[0];
sx q[0];
rz(0.71304148) q[0];
rz(2.9167922) q[1];
sx q[1];
rz(-0.59200042) q[1];
sx q[1];
rz(-1.0932659) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93738741) q[0];
sx q[0];
rz(-0.89417106) q[0];
sx q[0];
rz(2.1920374) q[0];
rz(2.1987783) q[2];
sx q[2];
rz(-0.33491376) q[2];
sx q[2];
rz(-0.57491344) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.374129) q[1];
sx q[1];
rz(-1.1157562) q[1];
sx q[1];
rz(-1.6947075) q[1];
rz(-1.1574) q[3];
sx q[3];
rz(-1.4634615) q[3];
sx q[3];
rz(-1.3232376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.0010058086) q[2];
sx q[2];
rz(-0.69585496) q[2];
sx q[2];
rz(0.33983964) q[2];
rz(-3.0991683) q[3];
sx q[3];
rz(-1.8570615) q[3];
sx q[3];
rz(-0.62780082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18702678) q[0];
sx q[0];
rz(-1.621959) q[0];
sx q[0];
rz(1.8931615) q[0];
rz(-0.78202248) q[1];
sx q[1];
rz(-0.93688688) q[1];
sx q[1];
rz(-0.72601906) q[1];
rz(1.9044884) q[2];
sx q[2];
rz(-1.2174229) q[2];
sx q[2];
rz(0.65897043) q[2];
rz(1.9966765) q[3];
sx q[3];
rz(-1.3419587) q[3];
sx q[3];
rz(1.6423196) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
