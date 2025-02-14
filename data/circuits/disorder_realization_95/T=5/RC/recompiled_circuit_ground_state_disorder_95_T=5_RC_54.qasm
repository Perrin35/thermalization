OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.65741444) q[0];
sx q[0];
rz(-1.8622458) q[0];
sx q[0];
rz(-0.6435464) q[0];
rz(0.14856385) q[1];
sx q[1];
rz(-2.4951976) q[1];
sx q[1];
rz(2.5656011) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51645378) q[0];
sx q[0];
rz(-1.0363434) q[0];
sx q[0];
rz(-2.778976) q[0];
rz(-1.2744489) q[2];
sx q[2];
rz(-1.7360032) q[2];
sx q[2];
rz(-2.8458386) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0514566) q[1];
sx q[1];
rz(-1.6701248) q[1];
sx q[1];
rz(-3.1121177) q[1];
rz(1.7891145) q[3];
sx q[3];
rz(-2.3872833) q[3];
sx q[3];
rz(1.393817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6286455) q[2];
sx q[2];
rz(-0.60428667) q[2];
sx q[2];
rz(-2.4224572) q[2];
rz(0.61270815) q[3];
sx q[3];
rz(-0.49379525) q[3];
sx q[3];
rz(-1.6814211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(1.7270108) q[0];
sx q[0];
rz(-2.5717323) q[0];
sx q[0];
rz(-1.3563096) q[0];
rz(2.5510229) q[1];
sx q[1];
rz(-1.5547724) q[1];
sx q[1];
rz(0.86033598) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3911678) q[0];
sx q[0];
rz(-2.2213687) q[0];
sx q[0];
rz(0.22919302) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27766547) q[2];
sx q[2];
rz(-0.57304731) q[2];
sx q[2];
rz(-0.78005698) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.4374784) q[1];
sx q[1];
rz(-1.8859004) q[1];
sx q[1];
rz(-1.9060288) q[1];
rz(-pi) q[2];
rz(1.4707056) q[3];
sx q[3];
rz(-2.2948569) q[3];
sx q[3];
rz(-2.049618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4071953) q[2];
sx q[2];
rz(-1.1043786) q[2];
sx q[2];
rz(2.8779136) q[2];
rz(-3.0401491) q[3];
sx q[3];
rz(-2.0439309) q[3];
sx q[3];
rz(1.4640079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2270671) q[0];
sx q[0];
rz(-2.3093746) q[0];
sx q[0];
rz(-1.9050003) q[0];
rz(-0.49691686) q[1];
sx q[1];
rz(-0.91941112) q[1];
sx q[1];
rz(0.20981728) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1001756) q[0];
sx q[0];
rz(-0.13209535) q[0];
sx q[0];
rz(2.9939674) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8489072) q[2];
sx q[2];
rz(-0.46574083) q[2];
sx q[2];
rz(-1.2276219) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4821149) q[1];
sx q[1];
rz(-2.6365197) q[1];
sx q[1];
rz(-1.494808) q[1];
x q[2];
rz(0.83720343) q[3];
sx q[3];
rz(-1.687369) q[3];
sx q[3];
rz(-1.0511412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0158374) q[2];
sx q[2];
rz(-2.0070984) q[2];
sx q[2];
rz(0.26975676) q[2];
rz(-0.45675993) q[3];
sx q[3];
rz(-0.41958198) q[3];
sx q[3];
rz(0.30341283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1064827) q[0];
sx q[0];
rz(-2.7484317) q[0];
sx q[0];
rz(-3.0495354) q[0];
rz(2.5606142) q[1];
sx q[1];
rz(-2.5216504) q[1];
sx q[1];
rz(-2.7041025) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70568181) q[0];
sx q[0];
rz(-1.4102524) q[0];
sx q[0];
rz(2.1401587) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3542489) q[2];
sx q[2];
rz(-1.4581973) q[2];
sx q[2];
rz(-1.6888022) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0713378) q[1];
sx q[1];
rz(-1.1769088) q[1];
sx q[1];
rz(2.8699234) q[1];
rz(-pi) q[2];
rz(-0.32467036) q[3];
sx q[3];
rz(-0.89327795) q[3];
sx q[3];
rz(3.0538525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.94365326) q[2];
sx q[2];
rz(-1.3760171) q[2];
sx q[2];
rz(0.1680689) q[2];
rz(-0.69593143) q[3];
sx q[3];
rz(-1.1974502) q[3];
sx q[3];
rz(-1.3958989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.593852) q[0];
sx q[0];
rz(-0.20741367) q[0];
sx q[0];
rz(-2.4399624) q[0];
rz(-0.54948366) q[1];
sx q[1];
rz(-2.0618849) q[1];
sx q[1];
rz(0.93542498) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9222636) q[0];
sx q[0];
rz(-0.50271243) q[0];
sx q[0];
rz(2.1070087) q[0];
rz(-pi) q[1];
rz(-2.7358908) q[2];
sx q[2];
rz(-3.0076111) q[2];
sx q[2];
rz(-1.9056428) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.11338698) q[1];
sx q[1];
rz(-2.1396239) q[1];
sx q[1];
rz(0.48888205) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0699877) q[3];
sx q[3];
rz(-1.6354899) q[3];
sx q[3];
rz(0.66401362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51297411) q[2];
sx q[2];
rz(-1.3209891) q[2];
sx q[2];
rz(-2.4690907) q[2];
rz(-2.3682112) q[3];
sx q[3];
rz(-2.0956109) q[3];
sx q[3];
rz(0.33236233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.8124939) q[0];
sx q[0];
rz(-2.2342873) q[0];
sx q[0];
rz(1.6727653) q[0];
rz(0.86917296) q[1];
sx q[1];
rz(-0.69907993) q[1];
sx q[1];
rz(2.8156143) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3968626) q[0];
sx q[0];
rz(-0.32856634) q[0];
sx q[0];
rz(2.7129125) q[0];
rz(-pi) q[1];
x q[1];
rz(2.432304) q[2];
sx q[2];
rz(-2.4457259) q[2];
sx q[2];
rz(2.8036731) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.80854844) q[1];
sx q[1];
rz(-2.1426629) q[1];
sx q[1];
rz(-0.99355917) q[1];
rz(1.4255089) q[3];
sx q[3];
rz(-0.43244468) q[3];
sx q[3];
rz(0.43859453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.088923067) q[2];
sx q[2];
rz(-1.576705) q[2];
sx q[2];
rz(2.3443473) q[2];
rz(3.0439324) q[3];
sx q[3];
rz(-0.74334136) q[3];
sx q[3];
rz(-1.9184448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.866211) q[0];
sx q[0];
rz(-2.1480063) q[0];
sx q[0];
rz(2.1589101) q[0];
rz(1.4879701) q[1];
sx q[1];
rz(-0.5653615) q[1];
sx q[1];
rz(2.8178094) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.569963) q[0];
sx q[0];
rz(-1.4870207) q[0];
sx q[0];
rz(2.8555968) q[0];
x q[1];
rz(2.1828142) q[2];
sx q[2];
rz(-1.806815) q[2];
sx q[2];
rz(2.4356206) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.52101719) q[1];
sx q[1];
rz(-1.6843101) q[1];
sx q[1];
rz(0.6599636) q[1];
rz(-pi) q[2];
rz(-1.9110095) q[3];
sx q[3];
rz(-2.0476626) q[3];
sx q[3];
rz(-2.8906156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4299778) q[2];
sx q[2];
rz(-0.51402503) q[2];
sx q[2];
rz(-0.035710486) q[2];
rz(-0.7401554) q[3];
sx q[3];
rz(-1.6868846) q[3];
sx q[3];
rz(1.9945701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78298727) q[0];
sx q[0];
rz(-2.7645223) q[0];
sx q[0];
rz(-0.20088917) q[0];
rz(-2.5783077) q[1];
sx q[1];
rz(-1.7726243) q[1];
sx q[1];
rz(-0.39931452) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4056643) q[0];
sx q[0];
rz(-0.20819323) q[0];
sx q[0];
rz(-2.3407779) q[0];
rz(-pi) q[1];
rz(0.041157965) q[2];
sx q[2];
rz(-2.0510457) q[2];
sx q[2];
rz(-0.89370773) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3914403) q[1];
sx q[1];
rz(-1.2160157) q[1];
sx q[1];
rz(-1.9998235) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7186728) q[3];
sx q[3];
rz(-2.3302571) q[3];
sx q[3];
rz(0.75046152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9575955) q[2];
sx q[2];
rz(-2.5338569) q[2];
sx q[2];
rz(0.11873928) q[2];
rz(0.27058288) q[3];
sx q[3];
rz(-0.99004531) q[3];
sx q[3];
rz(0.15429193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0305233) q[0];
sx q[0];
rz(-1.1708165) q[0];
sx q[0];
rz(0.4554553) q[0];
rz(-2.9083374) q[1];
sx q[1];
rz(-0.74359727) q[1];
sx q[1];
rz(0.0049678405) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99762395) q[0];
sx q[0];
rz(-1.652209) q[0];
sx q[0];
rz(-1.5366014) q[0];
rz(-pi) q[1];
rz(1.8894779) q[2];
sx q[2];
rz(-1.2010842) q[2];
sx q[2];
rz(-1.8869417) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4691866) q[1];
sx q[1];
rz(-0.71629706) q[1];
sx q[1];
rz(2.3138381) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6065323) q[3];
sx q[3];
rz(-0.4260582) q[3];
sx q[3];
rz(-0.96340513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3576971) q[2];
sx q[2];
rz(-1.8507439) q[2];
sx q[2];
rz(0.89007393) q[2];
rz(2.2642073) q[3];
sx q[3];
rz(-0.93534094) q[3];
sx q[3];
rz(-2.4963511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18168618) q[0];
sx q[0];
rz(-0.077597685) q[0];
sx q[0];
rz(2.087387) q[0];
rz(-2.8554845) q[1];
sx q[1];
rz(-2.095463) q[1];
sx q[1];
rz(-1.8792763) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6968203) q[0];
sx q[0];
rz(-0.69545924) q[0];
sx q[0];
rz(0.45375844) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20514078) q[2];
sx q[2];
rz(-1.8281422) q[2];
sx q[2];
rz(-0.3347476) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6773379) q[1];
sx q[1];
rz(-2.6644775) q[1];
sx q[1];
rz(2.7349763) q[1];
rz(-pi) q[2];
rz(-0.083857337) q[3];
sx q[3];
rz(-2.1241123) q[3];
sx q[3];
rz(0.38804873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.012764843) q[2];
sx q[2];
rz(-2.3598857) q[2];
sx q[2];
rz(-2.8774366) q[2];
rz(-2.9014897) q[3];
sx q[3];
rz(-1.0485579) q[3];
sx q[3];
rz(-0.31606328) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8964597) q[0];
sx q[0];
rz(-0.75440732) q[0];
sx q[0];
rz(-0.9040133) q[0];
rz(0.28884197) q[1];
sx q[1];
rz(-1.7542417) q[1];
sx q[1];
rz(-0.075686553) q[1];
rz(-0.31927989) q[2];
sx q[2];
rz(-1.3360595) q[2];
sx q[2];
rz(1.3315249) q[2];
rz(0.071500304) q[3];
sx q[3];
rz(-1.6115474) q[3];
sx q[3];
rz(-1.1781884) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
