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
rz(-2.3930385) q[0];
sx q[0];
rz(-1.8129803) q[0];
sx q[0];
rz(-0.42088977) q[0];
rz(-0.66840494) q[1];
sx q[1];
rz(-1.5937807) q[1];
sx q[1];
rz(-1.1512383) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65934718) q[0];
sx q[0];
rz(-1.3696241) q[0];
sx q[0];
rz(0.45659275) q[0];
rz(-pi) q[1];
rz(1.0082863) q[2];
sx q[2];
rz(-1.9676625) q[2];
sx q[2];
rz(-2.4728554) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.80025339) q[1];
sx q[1];
rz(-2.1113987) q[1];
sx q[1];
rz(-2.4990891) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9164768) q[3];
sx q[3];
rz(-0.68194032) q[3];
sx q[3];
rz(-0.038393858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0737754) q[2];
sx q[2];
rz(-0.99738085) q[2];
sx q[2];
rz(2.1201521) q[2];
rz(-1.3218309) q[3];
sx q[3];
rz(-2.6363711) q[3];
sx q[3];
rz(-2.5652313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8030871) q[0];
sx q[0];
rz(-1.0271238) q[0];
sx q[0];
rz(2.8402253) q[0];
rz(-2.001568) q[1];
sx q[1];
rz(-2.3224484) q[1];
sx q[1];
rz(-1.0637306) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6401705) q[0];
sx q[0];
rz(-0.29336624) q[0];
sx q[0];
rz(1.9984931) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3688886) q[2];
sx q[2];
rz(-0.81939745) q[2];
sx q[2];
rz(2.47987) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.3435522) q[1];
sx q[1];
rz(-2.4497089) q[1];
sx q[1];
rz(2.6380098) q[1];
x q[2];
rz(-3.0228457) q[3];
sx q[3];
rz(-2.1034965) q[3];
sx q[3];
rz(1.9783731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.86338824) q[2];
sx q[2];
rz(-0.74625838) q[2];
sx q[2];
rz(0.14872742) q[2];
rz(-1.9637828) q[3];
sx q[3];
rz(-1.9494467) q[3];
sx q[3];
rz(-1.8671794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3026368) q[0];
sx q[0];
rz(-1.9479072) q[0];
sx q[0];
rz(2.1394011) q[0];
rz(1.8110555) q[1];
sx q[1];
rz(-1.1794773) q[1];
sx q[1];
rz(2.5340714) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.73344) q[0];
sx q[0];
rz(-0.79793733) q[0];
sx q[0];
rz(-0.95592461) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.153468) q[2];
sx q[2];
rz(-0.4764423) q[2];
sx q[2];
rz(-0.63240766) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3866736) q[1];
sx q[1];
rz(-2.6133839) q[1];
sx q[1];
rz(0.22294238) q[1];
rz(-0.99276279) q[3];
sx q[3];
rz(-0.93621636) q[3];
sx q[3];
rz(1.6333333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9017631) q[2];
sx q[2];
rz(-0.87368691) q[2];
sx q[2];
rz(0.98881161) q[2];
rz(-1.6648939) q[3];
sx q[3];
rz(-1.8035382) q[3];
sx q[3];
rz(-0.57507676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.7779509) q[0];
sx q[0];
rz(-0.87894428) q[0];
sx q[0];
rz(-1.9212035) q[0];
rz(-1.3149423) q[1];
sx q[1];
rz(-1.5053791) q[1];
sx q[1];
rz(0.13994089) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8927887) q[0];
sx q[0];
rz(-0.99274764) q[0];
sx q[0];
rz(-1.0885574) q[0];
rz(-1.0331421) q[2];
sx q[2];
rz(-0.93574474) q[2];
sx q[2];
rz(2.6376574) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8257041) q[1];
sx q[1];
rz(-1.53535) q[1];
sx q[1];
rz(-2.0263544) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1871376) q[3];
sx q[3];
rz(-1.6373489) q[3];
sx q[3];
rz(0.47011791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4466897) q[2];
sx q[2];
rz(-1.0580772) q[2];
sx q[2];
rz(3.0042082) q[2];
rz(1.8047699) q[3];
sx q[3];
rz(-1.5788014) q[3];
sx q[3];
rz(-3.1409851) q[3];
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
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7430275) q[0];
sx q[0];
rz(-1.7337357) q[0];
sx q[0];
rz(-0.16125691) q[0];
rz(-2.297961) q[1];
sx q[1];
rz(-2.3944941) q[1];
sx q[1];
rz(-1.3074494) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6661412) q[0];
sx q[0];
rz(-1.3031811) q[0];
sx q[0];
rz(-2.3580736) q[0];
rz(-pi) q[1];
rz(1.2044123) q[2];
sx q[2];
rz(-1.7569524) q[2];
sx q[2];
rz(2.5103593) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54164179) q[1];
sx q[1];
rz(-0.58990084) q[1];
sx q[1];
rz(-1.4093424) q[1];
rz(-pi) q[2];
rz(0.67477711) q[3];
sx q[3];
rz(-1.0193079) q[3];
sx q[3];
rz(1.2434208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1696986) q[2];
sx q[2];
rz(-0.43115386) q[2];
sx q[2];
rz(2.2141854) q[2];
rz(-1.5291322) q[3];
sx q[3];
rz(-1.1011139) q[3];
sx q[3];
rz(-2.7839938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9277495) q[0];
sx q[0];
rz(-2.1723211) q[0];
sx q[0];
rz(1.7359605) q[0];
rz(1.4145781) q[1];
sx q[1];
rz(-2.3443293) q[1];
sx q[1];
rz(-0.79967156) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23223497) q[0];
sx q[0];
rz(-3.0785311) q[0];
sx q[0];
rz(-1.7574278) q[0];
x q[1];
rz(-0.39787103) q[2];
sx q[2];
rz(-2.8074773) q[2];
sx q[2];
rz(-1.0256922) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.97669125) q[1];
sx q[1];
rz(-1.8215239) q[1];
sx q[1];
rz(1.9216955) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.40414895) q[3];
sx q[3];
rz(-1.3938483) q[3];
sx q[3];
rz(-1.9623985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8863525) q[2];
sx q[2];
rz(-2.4031874) q[2];
sx q[2];
rz(-0.8026455) q[2];
rz(-2.8211527) q[3];
sx q[3];
rz(-1.8010062) q[3];
sx q[3];
rz(1.6360487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7771512) q[0];
sx q[0];
rz(-3.1138804) q[0];
sx q[0];
rz(0.030315422) q[0];
rz(1.8346571) q[1];
sx q[1];
rz(-2.0590643) q[1];
sx q[1];
rz(-0.62987769) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24014609) q[0];
sx q[0];
rz(-1.9909262) q[0];
sx q[0];
rz(1.3646558) q[0];
rz(-0.44158641) q[2];
sx q[2];
rz(-1.1418742) q[2];
sx q[2];
rz(-2.7457604) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3105433) q[1];
sx q[1];
rz(-1.8632952) q[1];
sx q[1];
rz(-1.7834375) q[1];
rz(-pi) q[2];
rz(1.685018) q[3];
sx q[3];
rz(-1.1291613) q[3];
sx q[3];
rz(0.87678443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1097172) q[2];
sx q[2];
rz(-0.66706359) q[2];
sx q[2];
rz(1.9943705) q[2];
rz(2.4540497) q[3];
sx q[3];
rz(-2.5496428) q[3];
sx q[3];
rz(-1.8182925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3654093) q[0];
sx q[0];
rz(-0.23615806) q[0];
sx q[0];
rz(-2.6614406) q[0];
rz(0.40223739) q[1];
sx q[1];
rz(-1.6419342) q[1];
sx q[1];
rz(0.38280907) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3671591) q[0];
sx q[0];
rz(-1.8978137) q[0];
sx q[0];
rz(0.22255442) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9391962) q[2];
sx q[2];
rz(-2.830626) q[2];
sx q[2];
rz(2.5045365) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9342707) q[1];
sx q[1];
rz(-1.8658966) q[1];
sx q[1];
rz(1.9701824) q[1];
rz(-1.6005101) q[3];
sx q[3];
rz(-1.7669356) q[3];
sx q[3];
rz(2.8063584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.39168921) q[2];
sx q[2];
rz(-0.73035556) q[2];
sx q[2];
rz(2.7833617) q[2];
rz(2.7273438) q[3];
sx q[3];
rz(-1.0284938) q[3];
sx q[3];
rz(-0.39345583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8527894) q[0];
sx q[0];
rz(-1.8599334) q[0];
sx q[0];
rz(2.7237256) q[0];
rz(1.4990384) q[1];
sx q[1];
rz(-2.7211029) q[1];
sx q[1];
rz(-2.5299759) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9754936) q[0];
sx q[0];
rz(-2.8320441) q[0];
sx q[0];
rz(-0.73348372) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5081399) q[2];
sx q[2];
rz(-1.776172) q[2];
sx q[2];
rz(-0.54772553) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.32595769) q[1];
sx q[1];
rz(-1.2919869) q[1];
sx q[1];
rz(-2.0626948) q[1];
rz(1.4889273) q[3];
sx q[3];
rz(-1.22504) q[3];
sx q[3];
rz(0.053178259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.186782) q[2];
sx q[2];
rz(-1.6839226) q[2];
sx q[2];
rz(2.7853454) q[2];
rz(0.48480836) q[3];
sx q[3];
rz(-1.7907413) q[3];
sx q[3];
rz(-3.1118605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.818882) q[0];
sx q[0];
rz(-1.1047381) q[0];
sx q[0];
rz(1.4622965) q[0];
rz(-2.7541584) q[1];
sx q[1];
rz(-2.2144364) q[1];
sx q[1];
rz(2.3364054) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.062267508) q[0];
sx q[0];
rz(-2.1342417) q[0];
sx q[0];
rz(-0.84765537) q[0];
rz(-pi) q[1];
rz(0.83909713) q[2];
sx q[2];
rz(-0.57575127) q[2];
sx q[2];
rz(-0.29345278) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.78811121) q[1];
sx q[1];
rz(-2.821021) q[1];
sx q[1];
rz(0.48133565) q[1];
rz(-1.8191387) q[3];
sx q[3];
rz(-1.9232009) q[3];
sx q[3];
rz(0.067841522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8055973) q[2];
sx q[2];
rz(-2.5184641) q[2];
sx q[2];
rz(1.222329) q[2];
rz(0.81021106) q[3];
sx q[3];
rz(-0.92456341) q[3];
sx q[3];
rz(-1.5626102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7280818) q[0];
sx q[0];
rz(-1.6376729) q[0];
sx q[0];
rz(-1.6117657) q[0];
rz(-1.275508) q[1];
sx q[1];
rz(-0.38560148) q[1];
sx q[1];
rz(1.7332981) q[1];
rz(-1.9553984) q[2];
sx q[2];
rz(-1.4438965) q[2];
sx q[2];
rz(0.039197103) q[2];
rz(0.90295193) q[3];
sx q[3];
rz(-1.9330238) q[3];
sx q[3];
rz(-0.15180363) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
