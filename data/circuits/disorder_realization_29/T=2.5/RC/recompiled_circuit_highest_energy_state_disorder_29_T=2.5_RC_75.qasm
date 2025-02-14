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
rz(1.2443378) q[0];
sx q[0];
rz(-0.79255784) q[0];
sx q[0];
rz(-0.25294024) q[0];
rz(2.3674372) q[1];
sx q[1];
rz(-2.7634662) q[1];
sx q[1];
rz(2.7546496) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8229503) q[0];
sx q[0];
rz(-1.3143263) q[0];
sx q[0];
rz(-1.9049113) q[0];
x q[1];
rz(-1.6574442) q[2];
sx q[2];
rz(-2.1468825) q[2];
sx q[2];
rz(2.224161) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7751037) q[1];
sx q[1];
rz(-1.4867884) q[1];
sx q[1];
rz(-1.4909527) q[1];
rz(-pi) q[2];
rz(1.0798595) q[3];
sx q[3];
rz(-1.0359952) q[3];
sx q[3];
rz(-2.8859069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.046254961) q[2];
sx q[2];
rz(-2.3825808) q[2];
sx q[2];
rz(-1.4026027) q[2];
rz(1.4143573) q[3];
sx q[3];
rz(-0.71310133) q[3];
sx q[3];
rz(-0.49899092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4614748) q[0];
sx q[0];
rz(-2.6533227) q[0];
sx q[0];
rz(-0.71037355) q[0];
rz(-1.0164227) q[1];
sx q[1];
rz(-1.2998394) q[1];
sx q[1];
rz(-0.76510915) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9020664) q[0];
sx q[0];
rz(-1.8886811) q[0];
sx q[0];
rz(-2.5796298) q[0];
rz(-pi) q[1];
rz(0.69074735) q[2];
sx q[2];
rz(-2.1155069) q[2];
sx q[2];
rz(-0.70554698) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.930508) q[1];
sx q[1];
rz(-0.74568664) q[1];
sx q[1];
rz(-1.0812283) q[1];
rz(1.9809057) q[3];
sx q[3];
rz(-2.7795459) q[3];
sx q[3];
rz(2.1459879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1171099) q[2];
sx q[2];
rz(-0.64579248) q[2];
sx q[2];
rz(-2.4369241) q[2];
rz(2.9674528) q[3];
sx q[3];
rz(-2.2691059) q[3];
sx q[3];
rz(-2.7599938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0891377) q[0];
sx q[0];
rz(-1.313504) q[0];
sx q[0];
rz(-2.254159) q[0];
rz(-1.2499836) q[1];
sx q[1];
rz(-1.2108112) q[1];
sx q[1];
rz(1.7807622) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.783238) q[0];
sx q[0];
rz(-1.1941507) q[0];
sx q[0];
rz(0.54623099) q[0];
x q[1];
rz(-3.0615882) q[2];
sx q[2];
rz(-1.7081877) q[2];
sx q[2];
rz(-0.72890711) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.20818921) q[1];
sx q[1];
rz(-1.6709242) q[1];
sx q[1];
rz(-0.80762819) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.181608) q[3];
sx q[3];
rz(-1.5312315) q[3];
sx q[3];
rz(-0.54007441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.36797324) q[2];
sx q[2];
rz(-2.7757288) q[2];
sx q[2];
rz(1.492929) q[2];
rz(2.6922373) q[3];
sx q[3];
rz(-1.2949233) q[3];
sx q[3];
rz(-1.2202643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0538977) q[0];
sx q[0];
rz(-0.59006417) q[0];
sx q[0];
rz(1.7328523) q[0];
rz(-2.159481) q[1];
sx q[1];
rz(-1.5928007) q[1];
sx q[1];
rz(1.4062448) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2236693) q[0];
sx q[0];
rz(-1.6817998) q[0];
sx q[0];
rz(1.3551606) q[0];
rz(-1.5179619) q[2];
sx q[2];
rz(-1.5312734) q[2];
sx q[2];
rz(-1.3827349) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.50735578) q[1];
sx q[1];
rz(-2.0869531) q[1];
sx q[1];
rz(0.28395758) q[1];
x q[2];
rz(0.26867433) q[3];
sx q[3];
rz(-2.2019535) q[3];
sx q[3];
rz(-2.1674066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1534319) q[2];
sx q[2];
rz(-2.1250696) q[2];
sx q[2];
rz(-0.15596685) q[2];
rz(-0.65034741) q[3];
sx q[3];
rz(-1.1740843) q[3];
sx q[3];
rz(-0.11421886) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0680189) q[0];
sx q[0];
rz(-1.8719712) q[0];
sx q[0];
rz(0.20075783) q[0];
rz(1.1772032) q[1];
sx q[1];
rz(-0.8526082) q[1];
sx q[1];
rz(1.9815365) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8667127) q[0];
sx q[0];
rz(-1.794941) q[0];
sx q[0];
rz(2.969541) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54520901) q[2];
sx q[2];
rz(-0.13449796) q[2];
sx q[2];
rz(-0.70869499) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9721891) q[1];
sx q[1];
rz(-2.6693925) q[1];
sx q[1];
rz(-0.71581033) q[1];
x q[2];
rz(-0.20108624) q[3];
sx q[3];
rz(-1.0498091) q[3];
sx q[3];
rz(0.21846427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3843627) q[2];
sx q[2];
rz(-2.195916) q[2];
sx q[2];
rz(-0.45822701) q[2];
rz(-0.630817) q[3];
sx q[3];
rz(-1.9061371) q[3];
sx q[3];
rz(2.6423776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.683627) q[0];
sx q[0];
rz(-0.85245913) q[0];
sx q[0];
rz(-0.0068579554) q[0];
rz(-2.9099756) q[1];
sx q[1];
rz(-1.3791142) q[1];
sx q[1];
rz(1.0622271) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0322965) q[0];
sx q[0];
rz(-1.8826906) q[0];
sx q[0];
rz(-1.7005141) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7825323) q[2];
sx q[2];
rz(-2.0604302) q[2];
sx q[2];
rz(-2.8608866) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0343218) q[1];
sx q[1];
rz(-0.91718972) q[1];
sx q[1];
rz(1.8922217) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3788578) q[3];
sx q[3];
rz(-0.71832359) q[3];
sx q[3];
rz(0.11107055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0007533) q[2];
sx q[2];
rz(-1.2522298) q[2];
sx q[2];
rz(-1.2678649) q[2];
rz(2.2800692) q[3];
sx q[3];
rz(-1.0637161) q[3];
sx q[3];
rz(1.4388194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
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
rz(-0.93290257) q[0];
sx q[0];
rz(-0.33354315) q[0];
sx q[0];
rz(-0.21555899) q[0];
rz(0.49628273) q[1];
sx q[1];
rz(-1.1880621) q[1];
sx q[1];
rz(-0.15942474) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8529786) q[0];
sx q[0];
rz(-0.63487494) q[0];
sx q[0];
rz(-1.575241) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47093289) q[2];
sx q[2];
rz(-1.3976946) q[2];
sx q[2];
rz(1.3927695) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7466399) q[1];
sx q[1];
rz(-1.8059019) q[1];
sx q[1];
rz(-0.04227738) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0219021) q[3];
sx q[3];
rz(-1.1682142) q[3];
sx q[3];
rz(0.91820133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.40519199) q[2];
sx q[2];
rz(-1.1242194) q[2];
sx q[2];
rz(0.51652017) q[2];
rz(-1.9308331) q[3];
sx q[3];
rz(-1.6885933) q[3];
sx q[3];
rz(-1.3794544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4844168) q[0];
sx q[0];
rz(-3.0653937) q[0];
sx q[0];
rz(-0.54022378) q[0];
rz(1.5243439) q[1];
sx q[1];
rz(-1.0018307) q[1];
sx q[1];
rz(1.2670955) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1138685) q[0];
sx q[0];
rz(-1.0719711) q[0];
sx q[0];
rz(-1.3214825) q[0];
rz(-pi) q[1];
rz(2.4300366) q[2];
sx q[2];
rz(-0.91469736) q[2];
sx q[2];
rz(-3.1392821) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9766751) q[1];
sx q[1];
rz(-1.7150214) q[1];
sx q[1];
rz(-0.37101908) q[1];
rz(-pi) q[2];
rz(-0.35532002) q[3];
sx q[3];
rz(-0.96002022) q[3];
sx q[3];
rz(-2.5988117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8871062) q[2];
sx q[2];
rz(-1.1595414) q[2];
sx q[2];
rz(-0.17733388) q[2];
rz(1.0698498) q[3];
sx q[3];
rz(-2.0900574) q[3];
sx q[3];
rz(2.9593318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7808481) q[0];
sx q[0];
rz(-1.1540664) q[0];
sx q[0];
rz(3.0408707) q[0];
rz(-1.1489457) q[1];
sx q[1];
rz(-1.9457685) q[1];
sx q[1];
rz(1.9675868) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7851104) q[0];
sx q[0];
rz(-2.504341) q[0];
sx q[0];
rz(-0.32109813) q[0];
rz(0.10020013) q[2];
sx q[2];
rz(-2.3728307) q[2];
sx q[2];
rz(-0.59610808) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9575038) q[1];
sx q[1];
rz(-3.0991249) q[1];
sx q[1];
rz(1.0951418) q[1];
rz(-pi) q[2];
rz(1.6948918) q[3];
sx q[3];
rz(-2.4467203) q[3];
sx q[3];
rz(-1.8852888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.164087) q[2];
sx q[2];
rz(-1.4848494) q[2];
sx q[2];
rz(2.738415) q[2];
rz(2.4781503) q[3];
sx q[3];
rz(-0.57938975) q[3];
sx q[3];
rz(0.0042393953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73096257) q[0];
sx q[0];
rz(-2.0616489) q[0];
sx q[0];
rz(3.0626815) q[0];
rz(2.2321189) q[1];
sx q[1];
rz(-1.0458922) q[1];
sx q[1];
rz(-0.39631072) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8412668) q[0];
sx q[0];
rz(-1.5551651) q[0];
sx q[0];
rz(-0.0090740694) q[0];
x q[1];
rz(-3.134507) q[2];
sx q[2];
rz(-0.79372294) q[2];
sx q[2];
rz(0.34863472) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9379314) q[1];
sx q[1];
rz(-2.0453369) q[1];
sx q[1];
rz(-1.306626) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.590606) q[3];
sx q[3];
rz(-2.1175623) q[3];
sx q[3];
rz(-1.3133698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4497455) q[2];
sx q[2];
rz(-2.2578466) q[2];
sx q[2];
rz(0.79908243) q[2];
rz(3.0633022) q[3];
sx q[3];
rz(-1.2543863) q[3];
sx q[3];
rz(-2.4962795) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66452022) q[0];
sx q[0];
rz(-1.399853) q[0];
sx q[0];
rz(0.54367263) q[0];
rz(-1.1935344) q[1];
sx q[1];
rz(-1.2763034) q[1];
sx q[1];
rz(-2.6313849) q[1];
rz(1.1134182) q[2];
sx q[2];
rz(-1.7241679) q[2];
sx q[2];
rz(0.67571251) q[2];
rz(2.067749) q[3];
sx q[3];
rz(-2.2324149) q[3];
sx q[3];
rz(0.80692337) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
