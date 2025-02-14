OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.34831369) q[0];
sx q[0];
rz(8.0170595) q[0];
sx q[0];
rz(9.4693139) q[0];
rz(2.0308004) q[1];
sx q[1];
rz(-1.9171311) q[1];
sx q[1];
rz(-2.3543624) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7101332) q[0];
sx q[0];
rz(-2.0636807) q[0];
sx q[0];
rz(-2.5566275) q[0];
rz(-pi) q[1];
rz(0.039723176) q[2];
sx q[2];
rz(-1.8958223) q[2];
sx q[2];
rz(-0.15234337) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8833784) q[1];
sx q[1];
rz(-1.6401263) q[1];
sx q[1];
rz(0.99154559) q[1];
x q[2];
rz(1.2177666) q[3];
sx q[3];
rz(-1.7952836) q[3];
sx q[3];
rz(0.29453555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.74701509) q[2];
sx q[2];
rz(-1.0801103) q[2];
sx q[2];
rz(0.75418312) q[2];
rz(1.7265823) q[3];
sx q[3];
rz(-0.49608803) q[3];
sx q[3];
rz(2.8126341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0542145) q[0];
sx q[0];
rz(-2.7890451) q[0];
sx q[0];
rz(-1.8074328) q[0];
rz(1.8502024) q[1];
sx q[1];
rz(-1.038237) q[1];
sx q[1];
rz(-0.87563595) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.310431) q[0];
sx q[0];
rz(-1.0410032) q[0];
sx q[0];
rz(-0.77052643) q[0];
rz(1.4562143) q[2];
sx q[2];
rz(-2.6422524) q[2];
sx q[2];
rz(-2.2407766) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3188556) q[1];
sx q[1];
rz(-1.8606288) q[1];
sx q[1];
rz(1.1136832) q[1];
x q[2];
rz(2.8505136) q[3];
sx q[3];
rz(-1.3716231) q[3];
sx q[3];
rz(1.8739623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4448173) q[2];
sx q[2];
rz(-2.4397218) q[2];
sx q[2];
rz(-0.87265054) q[2];
rz(-2.3526092) q[3];
sx q[3];
rz(-0.9011457) q[3];
sx q[3];
rz(-2.9130329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92662421) q[0];
sx q[0];
rz(-1.7890395) q[0];
sx q[0];
rz(-3.1335926) q[0];
rz(-2.3232715) q[1];
sx q[1];
rz(-0.74940959) q[1];
sx q[1];
rz(-1.9047033) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0950355) q[0];
sx q[0];
rz(-2.3716784) q[0];
sx q[0];
rz(-1.1004992) q[0];
rz(2.7929467) q[2];
sx q[2];
rz(-1.9816035) q[2];
sx q[2];
rz(-0.88627671) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9489386) q[1];
sx q[1];
rz(-2.1100419) q[1];
sx q[1];
rz(-1.9303028) q[1];
rz(0.45435702) q[3];
sx q[3];
rz(-2.8792692) q[3];
sx q[3];
rz(-1.9087877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7311953) q[2];
sx q[2];
rz(-0.13430139) q[2];
sx q[2];
rz(0.94089874) q[2];
rz(2.0375552) q[3];
sx q[3];
rz(-1.7887812) q[3];
sx q[3];
rz(-1.1436536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0333772) q[0];
sx q[0];
rz(-0.5431076) q[0];
sx q[0];
rz(2.9040842) q[0];
rz(-0.6595276) q[1];
sx q[1];
rz(-0.57412761) q[1];
sx q[1];
rz(3.1245756) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.983611) q[0];
sx q[0];
rz(-1.6189623) q[0];
sx q[0];
rz(-1.401813) q[0];
rz(-pi) q[1];
rz(0.41827664) q[2];
sx q[2];
rz(-2.661663) q[2];
sx q[2];
rz(-1.0470225) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.81777311) q[1];
sx q[1];
rz(-0.57460472) q[1];
sx q[1];
rz(-0.19293228) q[1];
rz(-pi) q[2];
rz(2.3364725) q[3];
sx q[3];
rz(-1.0110572) q[3];
sx q[3];
rz(0.92529682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.90618769) q[2];
sx q[2];
rz(-1.9402639) q[2];
sx q[2];
rz(2.2271633) q[2];
rz(-2.7794481) q[3];
sx q[3];
rz(-2.3320964) q[3];
sx q[3];
rz(-2.5563498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2528766) q[0];
sx q[0];
rz(-1.7676366) q[0];
sx q[0];
rz(-3.0294898) q[0];
rz(-1.0653227) q[1];
sx q[1];
rz(-1.0708829) q[1];
sx q[1];
rz(-1.0692474) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0651848) q[0];
sx q[0];
rz(-1.2994081) q[0];
sx q[0];
rz(-0.33581622) q[0];
x q[1];
rz(-2.6091796) q[2];
sx q[2];
rz(-1.9892009) q[2];
sx q[2];
rz(-2.6765228) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5204864) q[1];
sx q[1];
rz(-2.6992528) q[1];
sx q[1];
rz(-2.3250595) q[1];
x q[2];
rz(-1.3773243) q[3];
sx q[3];
rz(-2.5335823) q[3];
sx q[3];
rz(1.6093169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4273044) q[2];
sx q[2];
rz(-0.71791831) q[2];
sx q[2];
rz(0.53421268) q[2];
rz(-2.838375) q[3];
sx q[3];
rz(-1.5847619) q[3];
sx q[3];
rz(-1.5155972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7233906) q[0];
sx q[0];
rz(-1.8573107) q[0];
sx q[0];
rz(-2.2499625) q[0];
rz(1.3806237) q[1];
sx q[1];
rz(-1.9386407) q[1];
sx q[1];
rz(-0.92323971) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1311125) q[0];
sx q[0];
rz(-1.6690055) q[0];
sx q[0];
rz(1.6232071) q[0];
x q[1];
rz(-1.8452066) q[2];
sx q[2];
rz(-1.9447127) q[2];
sx q[2];
rz(-2.3309938) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6369317) q[1];
sx q[1];
rz(-1.8430222) q[1];
sx q[1];
rz(-2.9039761) q[1];
rz(-2.4999077) q[3];
sx q[3];
rz(-1.8998134) q[3];
sx q[3];
rz(-1.6365479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0455857) q[2];
sx q[2];
rz(-2.7092689) q[2];
sx q[2];
rz(-1.879479) q[2];
rz(1.9803842) q[3];
sx q[3];
rz(-1.5491756) q[3];
sx q[3];
rz(-1.8956634) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68814174) q[0];
sx q[0];
rz(-0.83204404) q[0];
sx q[0];
rz(1.6081109) q[0];
rz(-2.7104132) q[1];
sx q[1];
rz(-2.2240413) q[1];
sx q[1];
rz(1.4088438) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.726175) q[0];
sx q[0];
rz(-1.2841932) q[0];
sx q[0];
rz(2.4917401) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0529019) q[2];
sx q[2];
rz(-1.9541395) q[2];
sx q[2];
rz(2.2191522) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6366439) q[1];
sx q[1];
rz(-1.4145936) q[1];
sx q[1];
rz(1.4438932) q[1];
rz(-pi) q[2];
rz(-0.124229) q[3];
sx q[3];
rz(-1.0719187) q[3];
sx q[3];
rz(1.4829829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2048753) q[2];
sx q[2];
rz(-1.1254346) q[2];
sx q[2];
rz(-1.3196866) q[2];
rz(-3.1070869) q[3];
sx q[3];
rz(-1.6104108) q[3];
sx q[3];
rz(1.3694793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7337604) q[0];
sx q[0];
rz(-0.5540846) q[0];
sx q[0];
rz(1.9246509) q[0];
rz(2.2159684) q[1];
sx q[1];
rz(-1.3652912) q[1];
sx q[1];
rz(-1.2006753) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.417765) q[0];
sx q[0];
rz(-0.96207959) q[0];
sx q[0];
rz(-2.5136652) q[0];
rz(-pi) q[1];
rz(0.60007168) q[2];
sx q[2];
rz(-1.2286376) q[2];
sx q[2];
rz(-3.0759157) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4958805) q[1];
sx q[1];
rz(-0.79949841) q[1];
sx q[1];
rz(-1.6281566) q[1];
rz(-pi) q[2];
rz(1.993409) q[3];
sx q[3];
rz(-1.9070996) q[3];
sx q[3];
rz(2.5667532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.919148) q[2];
sx q[2];
rz(-1.5347975) q[2];
sx q[2];
rz(-2.0929125) q[2];
rz(0.18516304) q[3];
sx q[3];
rz(-1.984963) q[3];
sx q[3];
rz(3.0345501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1222526) q[0];
sx q[0];
rz(-2.4801319) q[0];
sx q[0];
rz(2.3308603) q[0];
rz(-0.73075378) q[1];
sx q[1];
rz(-2.0850875) q[1];
sx q[1];
rz(-0.96053851) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1584523) q[0];
sx q[0];
rz(-3.1380655) q[0];
sx q[0];
rz(-2.8615224) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4308295) q[2];
sx q[2];
rz(-1.1053876) q[2];
sx q[2];
rz(-0.39164603) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.62453485) q[1];
sx q[1];
rz(-1.0427999) q[1];
sx q[1];
rz(0.30450423) q[1];
rz(1.0455564) q[3];
sx q[3];
rz(-2.6572795) q[3];
sx q[3];
rz(-0.37176311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9591799) q[2];
sx q[2];
rz(-1.3591839) q[2];
sx q[2];
rz(-1.6020927) q[2];
rz(-1.0473853) q[3];
sx q[3];
rz(-1.3771219) q[3];
sx q[3];
rz(1.2654977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78028107) q[0];
sx q[0];
rz(-2.2762716) q[0];
sx q[0];
rz(-2.5541232) q[0];
rz(2.7777708) q[1];
sx q[1];
rz(-1.9963341) q[1];
sx q[1];
rz(-2.2344373) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1253352) q[0];
sx q[0];
rz(-1.2205219) q[0];
sx q[0];
rz(-1.6566234) q[0];
x q[1];
rz(-1.9586708) q[2];
sx q[2];
rz(-2.0128314) q[2];
sx q[2];
rz(1.3988053) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3876411) q[1];
sx q[1];
rz(-0.50261897) q[1];
sx q[1];
rz(2.6941264) q[1];
rz(-pi) q[2];
rz(1.7916405) q[3];
sx q[3];
rz(-2.7691602) q[3];
sx q[3];
rz(2.7999634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.63984799) q[2];
sx q[2];
rz(-3.0526243) q[2];
sx q[2];
rz(-2.7196344) q[2];
rz(0.0095327775) q[3];
sx q[3];
rz(-1.5671174) q[3];
sx q[3];
rz(2.6154521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6137153) q[0];
sx q[0];
rz(-1.8462702) q[0];
sx q[0];
rz(0.12311002) q[0];
rz(2.4404424) q[1];
sx q[1];
rz(-1.2260561) q[1];
sx q[1];
rz(-1.4969926) q[1];
rz(-0.43177615) q[2];
sx q[2];
rz(-0.74413055) q[2];
sx q[2];
rz(-1.6091138) q[2];
rz(-1.1604068) q[3];
sx q[3];
rz(-1.3153362) q[3];
sx q[3];
rz(1.5462331) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
