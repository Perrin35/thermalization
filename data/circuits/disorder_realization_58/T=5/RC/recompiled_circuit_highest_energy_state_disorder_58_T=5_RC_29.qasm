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
rz(0.27920029) q[0];
sx q[0];
rz(-1.2503799) q[0];
sx q[0];
rz(-0.00032902349) q[0];
rz(-2.012913) q[1];
sx q[1];
rz(-2.0560775) q[1];
sx q[1];
rz(0.20927277) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4930388) q[0];
sx q[0];
rz(-1.4161515) q[0];
sx q[0];
rz(-1.3376608) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5821788) q[2];
sx q[2];
rz(-1.7119679) q[2];
sx q[2];
rz(2.5641297) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.31262235) q[1];
sx q[1];
rz(-0.43401845) q[1];
sx q[1];
rz(2.390998) q[1];
x q[2];
rz(-2.9023323) q[3];
sx q[3];
rz(-1.3321028) q[3];
sx q[3];
rz(-2.8857114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2613232) q[2];
sx q[2];
rz(-1.3618733) q[2];
sx q[2];
rz(2.9051991) q[2];
rz(-2.7586625) q[3];
sx q[3];
rz(-1.087944) q[3];
sx q[3];
rz(1.6898164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7616538) q[0];
sx q[0];
rz(-0.28872696) q[0];
sx q[0];
rz(-2.5987103) q[0];
rz(0.44034964) q[1];
sx q[1];
rz(-2.0183225) q[1];
sx q[1];
rz(1.1276668) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(15/(8*pi)) q[0];
sx q[0];
rz(-2.4003599) q[0];
sx q[0];
rz(-2.6165589) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15717984) q[2];
sx q[2];
rz(-2.2008385) q[2];
sx q[2];
rz(-2.4147803) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.013036916) q[1];
sx q[1];
rz(-0.77289509) q[1];
sx q[1];
rz(-1.3580639) q[1];
rz(-pi) q[2];
rz(1.9502968) q[3];
sx q[3];
rz(-2.2436525) q[3];
sx q[3];
rz(1.4899474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8684964) q[2];
sx q[2];
rz(-2.345583) q[2];
sx q[2];
rz(-1.2924755) q[2];
rz(-0.83186197) q[3];
sx q[3];
rz(-1.4757194) q[3];
sx q[3];
rz(0.44124916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.454994) q[0];
sx q[0];
rz(-0.43231493) q[0];
sx q[0];
rz(1.5469714) q[0];
rz(3.1349685) q[1];
sx q[1];
rz(-2.6057656) q[1];
sx q[1];
rz(2.5033902) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6168882) q[0];
sx q[0];
rz(-2.0456438) q[0];
sx q[0];
rz(3.1010702) q[0];
x q[1];
rz(1.1630959) q[2];
sx q[2];
rz(-2.2981653) q[2];
sx q[2];
rz(3.0646366) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3824752) q[1];
sx q[1];
rz(-0.51847547) q[1];
sx q[1];
rz(2.3692959) q[1];
rz(-pi) q[2];
rz(2.9262337) q[3];
sx q[3];
rz(-1.0342068) q[3];
sx q[3];
rz(0.99765618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39614761) q[2];
sx q[2];
rz(-3.1148995) q[2];
sx q[2];
rz(2.2849582) q[2];
rz(2.4062697) q[3];
sx q[3];
rz(-1.7263128) q[3];
sx q[3];
rz(2.0935238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51744866) q[0];
sx q[0];
rz(-3.090591) q[0];
sx q[0];
rz(2.2080102) q[0];
rz(0.75992641) q[1];
sx q[1];
rz(-2.3257207) q[1];
sx q[1];
rz(1.5504206) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5869414) q[0];
sx q[0];
rz(-1.4822032) q[0];
sx q[0];
rz(0.2811801) q[0];
x q[1];
rz(-1.6408338) q[2];
sx q[2];
rz(-2.0067482) q[2];
sx q[2];
rz(0.34740123) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0821678) q[1];
sx q[1];
rz(-2.2520506) q[1];
sx q[1];
rz(-1.8685982) q[1];
x q[2];
rz(2.891734) q[3];
sx q[3];
rz(-1.7267531) q[3];
sx q[3];
rz(-0.74710959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7237225) q[2];
sx q[2];
rz(-1.7524065) q[2];
sx q[2];
rz(-2.1088481) q[2];
rz(-2.6241153) q[3];
sx q[3];
rz(-1.4314707) q[3];
sx q[3];
rz(1.3051858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1037769) q[0];
sx q[0];
rz(-3.0396437) q[0];
sx q[0];
rz(-1.284449) q[0];
rz(0.85353723) q[1];
sx q[1];
rz(-1.4505016) q[1];
sx q[1];
rz(-0.4761214) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42890756) q[0];
sx q[0];
rz(-1.0700912) q[0];
sx q[0];
rz(-0.34523532) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36557646) q[2];
sx q[2];
rz(-1.8387059) q[2];
sx q[2];
rz(0.6838201) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.13083982) q[1];
sx q[1];
rz(-0.71341842) q[1];
sx q[1];
rz(-2.6793007) q[1];
rz(-pi) q[2];
rz(0.83104844) q[3];
sx q[3];
rz(-1.5171192) q[3];
sx q[3];
rz(-0.37572467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7242754) q[2];
sx q[2];
rz(-0.83160916) q[2];
sx q[2];
rz(-1.6012491) q[2];
rz(3.0658412) q[3];
sx q[3];
rz(-1.3906393) q[3];
sx q[3];
rz(2.4116404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1002355) q[0];
sx q[0];
rz(-0.55158177) q[0];
sx q[0];
rz(1.3496189) q[0];
rz(1.2760108) q[1];
sx q[1];
rz(-0.7834692) q[1];
sx q[1];
rz(-1.0412201) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43364701) q[0];
sx q[0];
rz(-0.7793847) q[0];
sx q[0];
rz(-0.31275605) q[0];
x q[1];
rz(1.6068993) q[2];
sx q[2];
rz(-0.82335938) q[2];
sx q[2];
rz(-2.2359816) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.36633101) q[1];
sx q[1];
rz(-1.9420673) q[1];
sx q[1];
rz(-1.8410147) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8550664) q[3];
sx q[3];
rz(-2.3085924) q[3];
sx q[3];
rz(2.7018911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.87018806) q[2];
sx q[2];
rz(-1.1689593) q[2];
sx q[2];
rz(-2.1427587) q[2];
rz(2.6073604) q[3];
sx q[3];
rz(-1.9637008) q[3];
sx q[3];
rz(-1.6629201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6967195) q[0];
sx q[0];
rz(-0.98170009) q[0];
sx q[0];
rz(-1.5848507) q[0];
rz(-1.0818256) q[1];
sx q[1];
rz(-1.4993246) q[1];
sx q[1];
rz(-0.36344847) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2242607) q[0];
sx q[0];
rz(-1.5646569) q[0];
sx q[0];
rz(0.21256769) q[0];
rz(-pi) q[1];
rz(-2.8446782) q[2];
sx q[2];
rz(-1.468942) q[2];
sx q[2];
rz(-0.43019003) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.195943) q[1];
sx q[1];
rz(-0.87313014) q[1];
sx q[1];
rz(-0.5742258) q[1];
rz(-pi) q[2];
rz(2.0601252) q[3];
sx q[3];
rz(-2.8301468) q[3];
sx q[3];
rz(1.2490727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5536559) q[2];
sx q[2];
rz(-2.9755972) q[2];
sx q[2];
rz(1.0774277) q[2];
rz(1.59498) q[3];
sx q[3];
rz(-1.3225821) q[3];
sx q[3];
rz(2.4685629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3308291) q[0];
sx q[0];
rz(-1.6442693) q[0];
sx q[0];
rz(-0.28054917) q[0];
rz(2.5827017) q[1];
sx q[1];
rz(-2.2160896) q[1];
sx q[1];
rz(-1.9299318) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9292363) q[0];
sx q[0];
rz(-2.2693737) q[0];
sx q[0];
rz(-2.2792666) q[0];
rz(-pi) q[1];
rz(0.40139426) q[2];
sx q[2];
rz(-0.52782431) q[2];
sx q[2];
rz(-2.2032705) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.385541) q[1];
sx q[1];
rz(-1.5480528) q[1];
sx q[1];
rz(0.77994101) q[1];
x q[2];
rz(1.3727851) q[3];
sx q[3];
rz(-1.5757547) q[3];
sx q[3];
rz(-1.229242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.1739379) q[2];
sx q[2];
rz(-1.8123947) q[2];
sx q[2];
rz(2.3363414) q[2];
rz(-2.0179181) q[3];
sx q[3];
rz(-1.1403964) q[3];
sx q[3];
rz(2.8203188) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5241908) q[0];
sx q[0];
rz(-2.5179355) q[0];
sx q[0];
rz(-2.6362841) q[0];
rz(-2.4303923) q[1];
sx q[1];
rz(-2.2732747) q[1];
sx q[1];
rz(-2.0288846) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9049194) q[0];
sx q[0];
rz(-2.5917333) q[0];
sx q[0];
rz(2.5540094) q[0];
x q[1];
rz(2.4975487) q[2];
sx q[2];
rz(-0.13453787) q[2];
sx q[2];
rz(0.74468871) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1755274) q[1];
sx q[1];
rz(-1.6972145) q[1];
sx q[1];
rz(-2.011622) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6203513) q[3];
sx q[3];
rz(-0.55618868) q[3];
sx q[3];
rz(-2.9002528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.50735551) q[2];
sx q[2];
rz(-1.9017744) q[2];
sx q[2];
rz(-0.235454) q[2];
rz(-2.7638451) q[3];
sx q[3];
rz(-0.38769671) q[3];
sx q[3];
rz(-0.050067576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78558755) q[0];
sx q[0];
rz(-1.8536785) q[0];
sx q[0];
rz(-3.0872784) q[0];
rz(-2.3771225) q[1];
sx q[1];
rz(-0.46934325) q[1];
sx q[1];
rz(0.70942172) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55290907) q[0];
sx q[0];
rz(-1.8158661) q[0];
sx q[0];
rz(-0.10677307) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8807259) q[2];
sx q[2];
rz(-1.7013019) q[2];
sx q[2];
rz(2.5640287) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1349984) q[1];
sx q[1];
rz(-1.6853764) q[1];
sx q[1];
rz(-2.5496818) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4414431) q[3];
sx q[3];
rz(-1.23063) q[3];
sx q[3];
rz(2.8486855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8556122) q[2];
sx q[2];
rz(-1.6667655) q[2];
sx q[2];
rz(-0.093185514) q[2];
rz(0.5992254) q[3];
sx q[3];
rz(-0.56120509) q[3];
sx q[3];
rz(2.3563201) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7080606) q[0];
sx q[0];
rz(-0.96994937) q[0];
sx q[0];
rz(2.1920828) q[0];
rz(-3.0769729) q[1];
sx q[1];
rz(-0.043793543) q[1];
sx q[1];
rz(0.66088062) q[1];
rz(-1.4501889) q[2];
sx q[2];
rz(-0.55477886) q[2];
sx q[2];
rz(3.0045493) q[2];
rz(1.7568236) q[3];
sx q[3];
rz(-0.91384619) q[3];
sx q[3];
rz(-2.1887907) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
