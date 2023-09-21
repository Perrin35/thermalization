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
rz(2.9714669) q[0];
sx q[0];
rz(10.210769) q[0];
rz(0.6056447) q[1];
sx q[1];
rz(3.7925386) q[1];
sx q[1];
rz(8.7906919) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2274322) q[0];
sx q[0];
rz(-0.81613805) q[0];
sx q[0];
rz(3.1022275) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4889105) q[2];
sx q[2];
rz(-2.4856644) q[2];
sx q[2];
rz(-1.8413078) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7635599) q[1];
sx q[1];
rz(-2.0114007) q[1];
sx q[1];
rz(2.2256758) q[1];
x q[2];
rz(2.7414054) q[3];
sx q[3];
rz(-1.5923946) q[3];
sx q[3];
rz(-2.9633629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9156076) q[2];
sx q[2];
rz(-2.6245124) q[2];
sx q[2];
rz(-1.8784286) q[2];
rz(-1.2849215) q[3];
sx q[3];
rz(-1.6804755) q[3];
sx q[3];
rz(-2.8485956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9830575) q[0];
sx q[0];
rz(-1.2936658) q[0];
sx q[0];
rz(-0.43757004) q[0];
rz(-0.63105398) q[1];
sx q[1];
rz(-0.32866207) q[1];
sx q[1];
rz(-0.22110573) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9843922) q[0];
sx q[0];
rz(-1.385014) q[0];
sx q[0];
rz(-2.0322582) q[0];
rz(-1.2731304) q[2];
sx q[2];
rz(-2.7781099) q[2];
sx q[2];
rz(-3.1306981) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.21048966) q[1];
sx q[1];
rz(-2.0277129) q[1];
sx q[1];
rz(-1.4355852) q[1];
x q[2];
rz(-2.9651871) q[3];
sx q[3];
rz(-2.6325912) q[3];
sx q[3];
rz(-0.41364663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.21800403) q[2];
sx q[2];
rz(-1.4322586) q[2];
sx q[2];
rz(0.58829266) q[2];
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
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5730729) q[0];
sx q[0];
rz(-0.097671106) q[0];
sx q[0];
rz(-1.0004689) q[0];
rz(-2.4160642) q[1];
sx q[1];
rz(-2.0670481) q[1];
sx q[1];
rz(-2.3838938) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0453148) q[0];
sx q[0];
rz(-1.0490388) q[0];
sx q[0];
rz(1.9406712) q[0];
rz(-1.7311086) q[2];
sx q[2];
rz(-1.3659039) q[2];
sx q[2];
rz(-0.93970539) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.11848395) q[1];
sx q[1];
rz(-0.99220905) q[1];
sx q[1];
rz(-0.56225496) q[1];
rz(-pi) q[2];
rz(0.14962872) q[3];
sx q[3];
rz(-2.4976839) q[3];
sx q[3];
rz(3.0444037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2086601) q[2];
sx q[2];
rz(-2.9186086) q[2];
sx q[2];
rz(0.83646742) q[2];
rz(1.6992735) q[3];
sx q[3];
rz(-1.4740372) q[3];
sx q[3];
rz(-2.9377655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02012415) q[0];
sx q[0];
rz(-0.080105372) q[0];
sx q[0];
rz(0.51112038) q[0];
rz(-0.077443667) q[1];
sx q[1];
rz(-0.42235342) q[1];
sx q[1];
rz(1.8130594) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0275832) q[0];
sx q[0];
rz(-2.8767577) q[0];
sx q[0];
rz(-0.94075216) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6287597) q[2];
sx q[2];
rz(-2.3094258) q[2];
sx q[2];
rz(1.8682478) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36973876) q[1];
sx q[1];
rz(-2.8747897) q[1];
sx q[1];
rz(-1.5727444) q[1];
rz(2.02416) q[3];
sx q[3];
rz(-0.24586596) q[3];
sx q[3];
rz(2.9256431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.093734309) q[2];
sx q[2];
rz(-1.2183943) q[2];
sx q[2];
rz(2.034534) q[2];
rz(2.6691061) q[3];
sx q[3];
rz(-1.5269591) q[3];
sx q[3];
rz(-0.85737491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1012786) q[0];
sx q[0];
rz(-2.2450228) q[0];
sx q[0];
rz(0.3381981) q[0];
rz(-1.2942554) q[1];
sx q[1];
rz(-1.5816403) q[1];
sx q[1];
rz(0.50450528) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11855928) q[0];
sx q[0];
rz(-2.1001864) q[0];
sx q[0];
rz(-2.8708007) q[0];
x q[1];
rz(1.9485103) q[2];
sx q[2];
rz(-1.1337122) q[2];
sx q[2];
rz(1.2010241) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.63197631) q[1];
sx q[1];
rz(-0.39034931) q[1];
sx q[1];
rz(2.2401287) q[1];
rz(-0.91315956) q[3];
sx q[3];
rz(-1.0201766) q[3];
sx q[3];
rz(-0.0084358128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.111104) q[2];
sx q[2];
rz(-2.3031074) q[2];
sx q[2];
rz(1.6476691) q[2];
rz(-1.4533639) q[3];
sx q[3];
rz(-2.2036392) q[3];
sx q[3];
rz(1.3822314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9313653) q[0];
sx q[0];
rz(-2.2670822) q[0];
sx q[0];
rz(-2.0507623) q[0];
rz(2.6018654) q[1];
sx q[1];
rz(-1.153774) q[1];
sx q[1];
rz(2.9528023) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8180346) q[0];
sx q[0];
rz(-1.8206882) q[0];
sx q[0];
rz(0.37139335) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9300869) q[2];
sx q[2];
rz(-0.28887666) q[2];
sx q[2];
rz(-1.2668244) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.530045) q[1];
sx q[1];
rz(-1.6366742) q[1];
sx q[1];
rz(-0.24286119) q[1];
x q[2];
rz(-0.24174989) q[3];
sx q[3];
rz(-1.3931735) q[3];
sx q[3];
rz(1.0279946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9699036) q[2];
sx q[2];
rz(-1.0287372) q[2];
sx q[2];
rz(1.8264654) q[2];
rz(1.5054437) q[3];
sx q[3];
rz(-2.7841778) q[3];
sx q[3];
rz(1.858254) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.090102) q[0];
sx q[0];
rz(-2.3110456) q[0];
sx q[0];
rz(2.8175957) q[0];
rz(1.8404768) q[1];
sx q[1];
rz(-2.3062861) q[1];
sx q[1];
rz(-1.7623998) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1349072) q[0];
sx q[0];
rz(-0.25776699) q[0];
sx q[0];
rz(-1.8968614) q[0];
x q[1];
rz(1.1629282) q[2];
sx q[2];
rz(-1.1424354) q[2];
sx q[2];
rz(0.63901627) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6143601) q[1];
sx q[1];
rz(-0.34253201) q[1];
sx q[1];
rz(-0.50369461) q[1];
x q[2];
rz(0.58052766) q[3];
sx q[3];
rz(-2.1248098) q[3];
sx q[3];
rz(-1.1298657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1088915) q[2];
sx q[2];
rz(-1.4543433) q[2];
sx q[2];
rz(0.94318715) q[2];
rz(-0.33106783) q[3];
sx q[3];
rz(-1.3676684) q[3];
sx q[3];
rz(0.29512063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11809764) q[0];
sx q[0];
rz(-0.58693111) q[0];
sx q[0];
rz(-2.3773637) q[0];
rz(-3.0006192) q[1];
sx q[1];
rz(-0.39607513) q[1];
sx q[1];
rz(-1.8483298) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4828669) q[0];
sx q[0];
rz(-1.3728766) q[0];
sx q[0];
rz(1.1307964) q[0];
rz(-0.60111945) q[2];
sx q[2];
rz(-1.4919315) q[2];
sx q[2];
rz(2.2476946) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5544719) q[1];
sx q[1];
rz(-1.1677824) q[1];
sx q[1];
rz(0.15092571) q[1];
rz(-2.7420298) q[3];
sx q[3];
rz(-1.385653) q[3];
sx q[3];
rz(-0.012133908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3747037) q[2];
sx q[2];
rz(-2.1024487) q[2];
sx q[2];
rz(-2.5602706) q[2];
rz(-2.2733722) q[3];
sx q[3];
rz(-2.7243082) q[3];
sx q[3];
rz(-1.3114312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.593489) q[0];
sx q[0];
rz(-0.5797525) q[0];
sx q[0];
rz(-1.2506437) q[0];
rz(0.99682322) q[1];
sx q[1];
rz(-0.88368982) q[1];
sx q[1];
rz(1.2876127) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2244789) q[0];
sx q[0];
rz(-2.0600442) q[0];
sx q[0];
rz(3.1130303) q[0];
x q[1];
rz(1.7436696) q[2];
sx q[2];
rz(-1.0908974) q[2];
sx q[2];
rz(-2.929504) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1907562) q[1];
sx q[1];
rz(-1.3637929) q[1];
sx q[1];
rz(0.91916577) q[1];
rz(-0.82448126) q[3];
sx q[3];
rz(-1.5319676) q[3];
sx q[3];
rz(-0.83394921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8210956) q[2];
sx q[2];
rz(-2.7536776) q[2];
sx q[2];
rz(-1.8851177) q[2];
rz(-2.9750032) q[3];
sx q[3];
rz(-1.5766671) q[3];
sx q[3];
rz(2.0690209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(-2.4644641) q[1];
sx q[1];
rz(-0.36718711) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31736483) q[0];
sx q[0];
rz(-1.6019078) q[0];
sx q[0];
rz(3.0935862) q[0];
rz(-pi) q[1];
rz(-1.2904097) q[2];
sx q[2];
rz(-0.9624316) q[2];
sx q[2];
rz(-0.3384564) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.90875188) q[1];
sx q[1];
rz(-2.4526261) q[1];
sx q[1];
rz(2.8244551) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.48093421) q[3];
sx q[3];
rz(-2.1554865) q[3];
sx q[3];
rz(1.9422873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.24511589) q[2];
sx q[2];
rz(-0.8090691) q[2];
sx q[2];
rz(1.0661351) q[2];
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
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29522482) q[0];
sx q[0];
rz(-1.5867148) q[0];
sx q[0];
rz(-1.5750194) q[0];
rz(-0.29905839) q[1];
sx q[1];
rz(-0.60332861) q[1];
sx q[1];
rz(0.52287846) q[1];
rz(-3.0565312) q[2];
sx q[2];
rz(-0.75123514) q[2];
sx q[2];
rz(-3.0820465) q[2];
rz(1.7442262) q[3];
sx q[3];
rz(-2.4122824) q[3];
sx q[3];
rz(-1.1576049) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];