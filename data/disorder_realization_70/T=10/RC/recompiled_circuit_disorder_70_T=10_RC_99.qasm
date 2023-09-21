OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.45286173) q[0];
sx q[0];
rz(-0.17012574) q[0];
sx q[0];
rz(-0.78599077) q[0];
rz(-2.535948) q[1];
sx q[1];
rz(-0.65094596) q[1];
sx q[1];
rz(-2.5075066) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4579826) q[0];
sx q[0];
rz(-1.542122) q[0];
sx q[0];
rz(2.3258414) q[0];
rz(-1.4889105) q[2];
sx q[2];
rz(-0.6559283) q[2];
sx q[2];
rz(1.3002849) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.69971985) q[1];
sx q[1];
rz(-2.3708214) q[1];
sx q[1];
rz(-0.91200478) q[1];
rz(-0.4001873) q[3];
sx q[3];
rz(-1.5923946) q[3];
sx q[3];
rz(0.17822972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9156076) q[2];
sx q[2];
rz(-2.6245124) q[2];
sx q[2];
rz(-1.263164) q[2];
rz(-1.2849215) q[3];
sx q[3];
rz(-1.6804755) q[3];
sx q[3];
rz(-2.8485956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9830575) q[0];
sx q[0];
rz(-1.2936658) q[0];
sx q[0];
rz(0.43757004) q[0];
rz(-2.5105387) q[1];
sx q[1];
rz(-0.32866207) q[1];
sx q[1];
rz(-2.9204869) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32199931) q[0];
sx q[0];
rz(-2.0237192) q[0];
sx q[0];
rz(2.9346912) q[0];
rz(-0.11110335) q[2];
sx q[2];
rz(-1.2239893) q[2];
sx q[2];
rz(2.8135516) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7213388) q[1];
sx q[1];
rz(-1.6920648) q[1];
sx q[1];
rz(2.6810357) q[1];
rz(-pi) q[2];
rz(2.6392194) q[3];
sx q[3];
rz(-1.4851735) q[3];
sx q[3];
rz(1.3115713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.21800403) q[2];
sx q[2];
rz(-1.709334) q[2];
sx q[2];
rz(2.5533) q[2];
rz(-2.6925987) q[3];
sx q[3];
rz(-0.42789999) q[3];
sx q[3];
rz(-2.0935521) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56851971) q[0];
sx q[0];
rz(-3.0439215) q[0];
sx q[0];
rz(-1.0004689) q[0];
rz(2.4160642) q[1];
sx q[1];
rz(-2.0670481) q[1];
sx q[1];
rz(-0.75769889) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.706447) q[0];
sx q[0];
rz(-2.5120814) q[0];
sx q[0];
rz(-0.56133095) q[0];
x q[1];
rz(-2.9341142) q[2];
sx q[2];
rz(-1.7277272) q[2];
sx q[2];
rz(-0.66397882) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.12049) q[1];
sx q[1];
rz(-2.0334525) q[1];
sx q[1];
rz(-0.91336577) q[1];
rz(-1.682231) q[3];
sx q[3];
rz(-0.93524869) q[3];
sx q[3];
rz(-0.089126822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9329325) q[2];
sx q[2];
rz(-0.22298403) q[2];
sx q[2];
rz(2.3051252) q[2];
rz(1.4423192) q[3];
sx q[3];
rz(-1.4740372) q[3];
sx q[3];
rz(2.9377655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02012415) q[0];
sx q[0];
rz(-3.0614873) q[0];
sx q[0];
rz(-0.51112038) q[0];
rz(3.064149) q[1];
sx q[1];
rz(-0.42235342) q[1];
sx q[1];
rz(-1.3285332) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1140095) q[0];
sx q[0];
rz(-2.8767577) q[0];
sx q[0];
rz(0.94075216) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51283299) q[2];
sx q[2];
rz(-2.3094258) q[2];
sx q[2];
rz(1.8682478) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9424142) q[1];
sx q[1];
rz(-1.5702827) q[1];
sx q[1];
rz(-1.3039939) q[1];
rz(-pi) q[2];
rz(-3.0321211) q[3];
sx q[3];
rz(-1.3502035) q[3];
sx q[3];
rz(0.68144875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.093734309) q[2];
sx q[2];
rz(-1.9231984) q[2];
sx q[2];
rz(-1.1070586) q[2];
rz(0.47248653) q[3];
sx q[3];
rz(-1.6146336) q[3];
sx q[3];
rz(-0.85737491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-3.1012786) q[0];
sx q[0];
rz(-0.89656985) q[0];
sx q[0];
rz(-0.3381981) q[0];
rz(1.8473373) q[1];
sx q[1];
rz(-1.5816403) q[1];
sx q[1];
rz(0.50450528) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11855928) q[0];
sx q[0];
rz(-1.0414062) q[0];
sx q[0];
rz(-2.8708007) q[0];
x q[1];
rz(1.9485103) q[2];
sx q[2];
rz(-2.0078805) q[2];
sx q[2];
rz(-1.2010241) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5096163) q[1];
sx q[1];
rz(-0.39034931) q[1];
sx q[1];
rz(0.901464) q[1];
rz(-pi) q[2];
rz(-0.91315956) q[3];
sx q[3];
rz(-1.0201766) q[3];
sx q[3];
rz(3.1331568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0304886) q[2];
sx q[2];
rz(-0.8384853) q[2];
sx q[2];
rz(1.6476691) q[2];
rz(-1.4533639) q[3];
sx q[3];
rz(-2.2036392) q[3];
sx q[3];
rz(-1.7593613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(2.9313653) q[0];
sx q[0];
rz(-2.2670822) q[0];
sx q[0];
rz(1.0908303) q[0];
rz(-0.53972721) q[1];
sx q[1];
rz(-1.153774) q[1];
sx q[1];
rz(-0.18879034) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81297368) q[0];
sx q[0];
rz(-2.6972174) q[0];
sx q[0];
rz(0.61291738) q[0];
rz(-1.2994453) q[2];
sx q[2];
rz(-1.6711298) q[2];
sx q[2];
rz(0.041610418) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.61154762) q[1];
sx q[1];
rz(-1.6366742) q[1];
sx q[1];
rz(-0.24286119) q[1];
rz(-pi) q[2];
rz(-0.64340274) q[3];
sx q[3];
rz(-2.8426369) q[3];
sx q[3];
rz(1.1645731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.171689) q[2];
sx q[2];
rz(-2.1128555) q[2];
sx q[2];
rz(1.8264654) q[2];
rz(-1.5054437) q[3];
sx q[3];
rz(-2.7841778) q[3];
sx q[3];
rz(1.2833387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0514907) q[0];
sx q[0];
rz(-0.83054709) q[0];
sx q[0];
rz(-2.8175957) q[0];
rz(1.8404768) q[1];
sx q[1];
rz(-2.3062861) q[1];
sx q[1];
rz(1.3791929) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8934879) q[0];
sx q[0];
rz(-1.4890492) q[0];
sx q[0];
rz(1.8155314) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6799913) q[2];
sx q[2];
rz(-1.9399376) q[2];
sx q[2];
rz(-0.75418562) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5768347) q[1];
sx q[1];
rz(-1.733629) q[1];
sx q[1];
rz(-2.8388883) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2961388) q[3];
sx q[3];
rz(-0.77973706) q[3];
sx q[3];
rz(-0.23507915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1088915) q[2];
sx q[2];
rz(-1.6872493) q[2];
sx q[2];
rz(2.1984055) q[2];
rz(2.8105248) q[3];
sx q[3];
rz(-1.3676684) q[3];
sx q[3];
rz(-2.846472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.023495) q[0];
sx q[0];
rz(-2.5546615) q[0];
sx q[0];
rz(-0.76422894) q[0];
rz(0.14097342) q[1];
sx q[1];
rz(-0.39607513) q[1];
sx q[1];
rz(1.2932628) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9613567) q[0];
sx q[0];
rz(-2.0016252) q[0];
sx q[0];
rz(0.21813099) q[0];
rz(-pi) q[1];
rz(1.4752611) q[2];
sx q[2];
rz(-0.97180688) q[2];
sx q[2];
rz(2.4107188) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0656933) q[1];
sx q[1];
rz(-1.7095487) q[1];
sx q[1];
rz(1.1636415) q[1];
x q[2];
rz(-0.39956283) q[3];
sx q[3];
rz(-1.7559397) q[3];
sx q[3];
rz(-0.012133908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3747037) q[2];
sx q[2];
rz(-2.1024487) q[2];
sx q[2];
rz(0.58132201) q[2];
rz(0.86822048) q[3];
sx q[3];
rz(-2.7243082) q[3];
sx q[3];
rz(1.8301615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.593489) q[0];
sx q[0];
rz(-2.5618401) q[0];
sx q[0];
rz(-1.8909489) q[0];
rz(-2.1447694) q[1];
sx q[1];
rz(-0.88368982) q[1];
sx q[1];
rz(1.2876127) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33289136) q[0];
sx q[0];
rz(-1.5455855) q[0];
sx q[0];
rz(2.0602134) q[0];
rz(-pi) q[1];
rz(-2.655517) q[2];
sx q[2];
rz(-1.7239778) q[2];
sx q[2];
rz(1.7024405) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3661256) q[1];
sx q[1];
rz(-2.2062416) q[1];
sx q[1];
rz(0.25823621) q[1];
rz(-pi) q[2];
rz(3.0887293) q[3];
sx q[3];
rz(-2.3164146) q[3];
sx q[3];
rz(2.4406274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8210956) q[2];
sx q[2];
rz(-0.38791502) q[2];
sx q[2];
rz(-1.8851177) q[2];
rz(-0.16658941) q[3];
sx q[3];
rz(-1.5766671) q[3];
sx q[3];
rz(-2.0690209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2607516) q[0];
sx q[0];
rz(-2.3374225) q[0];
sx q[0];
rz(0.32521954) q[0];
rz(2.0064158) q[1];
sx q[1];
rz(-0.67712855) q[1];
sx q[1];
rz(-2.7744055) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8866667) q[0];
sx q[0];
rz(-1.6187795) q[0];
sx q[0];
rz(-1.6019437) q[0];
rz(-2.514421) q[2];
sx q[2];
rz(-1.7998724) q[2];
sx q[2];
rz(1.0692182) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.91017427) q[1];
sx q[1];
rz(-1.3712198) q[1];
sx q[1];
rz(2.477596) q[1];
x q[2];
rz(-0.96079798) q[3];
sx q[3];
rz(-2.4028117) q[3];
sx q[3];
rz(1.184954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8964768) q[2];
sx q[2];
rz(-0.8090691) q[2];
sx q[2];
rz(-2.0754576) q[2];
rz(-3.0623479) q[3];
sx q[3];
rz(-2.5388122) q[3];
sx q[3];
rz(-1.4762896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29522482) q[0];
sx q[0];
rz(-1.5548779) q[0];
sx q[0];
rz(1.5665733) q[0];
rz(0.29905839) q[1];
sx q[1];
rz(-2.538264) q[1];
sx q[1];
rz(-2.6187142) q[1];
rz(3.0565312) q[2];
sx q[2];
rz(-2.3903575) q[2];
sx q[2];
rz(0.059546197) q[2];
rz(-0.15300898) q[3];
sx q[3];
rz(-2.2867793) q[3];
sx q[3];
rz(2.2147562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];