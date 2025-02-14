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
rz(0.74855411) q[0];
sx q[0];
rz(-1.3286123) q[0];
sx q[0];
rz(-2.7207029) q[0];
rz(2.4731877) q[1];
sx q[1];
rz(-1.547812) q[1];
sx q[1];
rz(1.1512383) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0092888) q[0];
sx q[0];
rz(-2.0175066) q[0];
sx q[0];
rz(1.7942091) q[0];
rz(0.45999476) q[2];
sx q[2];
rz(-1.056571) q[2];
sx q[2];
rz(-2.0005039) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0033773) q[1];
sx q[1];
rz(-1.0312407) q[1];
sx q[1];
rz(-0.92745499) q[1];
x q[2];
rz(0.26845308) q[3];
sx q[3];
rz(-2.2055948) q[3];
sx q[3];
rz(2.6688547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0737754) q[2];
sx q[2];
rz(-0.99738085) q[2];
sx q[2];
rz(1.0214405) q[2];
rz(1.8197618) q[3];
sx q[3];
rz(-0.50522155) q[3];
sx q[3];
rz(-0.57636133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
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
rz(-0.3385056) q[0];
sx q[0];
rz(-1.0271238) q[0];
sx q[0];
rz(-2.8402253) q[0];
rz(1.1400247) q[1];
sx q[1];
rz(-2.3224484) q[1];
sx q[1];
rz(-1.0637306) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0845818) q[0];
sx q[0];
rz(-1.8370596) q[0];
sx q[0];
rz(-0.12464704) q[0];
rz(-pi) q[1];
rz(1.3688886) q[2];
sx q[2];
rz(-2.3221952) q[2];
sx q[2];
rz(-2.47987) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.96460198) q[1];
sx q[1];
rz(-2.1637205) q[1];
sx q[1];
rz(-1.1904535) q[1];
x q[2];
rz(1.3724907) q[3];
sx q[3];
rz(-0.54452678) q[3];
sx q[3];
rz(0.93248065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.86338824) q[2];
sx q[2];
rz(-2.3953343) q[2];
sx q[2];
rz(-0.14872742) q[2];
rz(-1.1778098) q[3];
sx q[3];
rz(-1.1921459) q[3];
sx q[3];
rz(-1.8671794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.3026368) q[0];
sx q[0];
rz(-1.9479072) q[0];
sx q[0];
rz(1.0021915) q[0];
rz(-1.8110555) q[1];
sx q[1];
rz(-1.9621153) q[1];
sx q[1];
rz(-0.60752121) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70456767) q[0];
sx q[0];
rz(-1.9965197) q[0];
sx q[0];
rz(-0.87354891) q[0];
rz(0.98812466) q[2];
sx q[2];
rz(-2.6651504) q[2];
sx q[2];
rz(0.63240766) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5190815) q[1];
sx q[1];
rz(-1.6824598) q[1];
sx q[1];
rz(-2.6242248) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5030977) q[3];
sx q[3];
rz(-0.83052626) q[3];
sx q[3];
rz(2.3414224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2398296) q[2];
sx q[2];
rz(-0.87368691) q[2];
sx q[2];
rz(-0.98881161) q[2];
rz(-1.4766988) q[3];
sx q[3];
rz(-1.3380545) q[3];
sx q[3];
rz(-0.57507676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3636417) q[0];
sx q[0];
rz(-2.2626484) q[0];
sx q[0];
rz(-1.2203891) q[0];
rz(-1.3149423) q[1];
sx q[1];
rz(-1.6362135) q[1];
sx q[1];
rz(3.0016518) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24880399) q[0];
sx q[0];
rz(-0.99274764) q[0];
sx q[0];
rz(2.0530353) q[0];
x q[1];
rz(-2.432517) q[2];
sx q[2];
rz(-1.1458414) q[2];
sx q[2];
rz(2.4147025) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.23754696) q[1];
sx q[1];
rz(-2.0260467) q[1];
sx q[1];
rz(-3.1021248) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3945822) q[3];
sx q[3];
rz(-0.38910633) q[3];
sx q[3];
rz(-1.8776304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6949029) q[2];
sx q[2];
rz(-2.0835154) q[2];
sx q[2];
rz(0.13738446) q[2];
rz(-1.8047699) q[3];
sx q[3];
rz(-1.5788014) q[3];
sx q[3];
rz(3.1409851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7430275) q[0];
sx q[0];
rz(-1.7337357) q[0];
sx q[0];
rz(-2.9803357) q[0];
rz(2.297961) q[1];
sx q[1];
rz(-0.7470986) q[1];
sx q[1];
rz(1.8341433) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6661412) q[0];
sx q[0];
rz(-1.8384116) q[0];
sx q[0];
rz(-2.3580736) q[0];
x q[1];
rz(1.9371804) q[2];
sx q[2];
rz(-1.7569524) q[2];
sx q[2];
rz(-2.5103593) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89462639) q[1];
sx q[1];
rz(-1.481253) q[1];
sx q[1];
rz(-0.98693165) q[1];
rz(-pi) q[2];
rz(-2.3639115) q[3];
sx q[3];
rz(-2.2984004) q[3];
sx q[3];
rz(2.889159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1696986) q[2];
sx q[2];
rz(-0.43115386) q[2];
sx q[2];
rz(2.2141854) q[2];
rz(1.5291322) q[3];
sx q[3];
rz(-2.0404787) q[3];
sx q[3];
rz(-2.7839938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2138432) q[0];
sx q[0];
rz(-2.1723211) q[0];
sx q[0];
rz(1.4056322) q[0];
rz(-1.7270145) q[1];
sx q[1];
rz(-2.3443293) q[1];
sx q[1];
rz(2.3419211) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9093577) q[0];
sx q[0];
rz(-3.0785311) q[0];
sx q[0];
rz(-1.3841649) q[0];
rz(-pi) q[1];
rz(1.437101) q[2];
sx q[2];
rz(-1.8779199) q[2];
sx q[2];
rz(-1.6973059) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.97669125) q[1];
sx q[1];
rz(-1.8215239) q[1];
sx q[1];
rz(-1.2198971) q[1];
rz(-pi) q[2];
x q[2];
rz(1.762885) q[3];
sx q[3];
rz(-1.1733161) q[3];
sx q[3];
rz(2.8251338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.25524011) q[2];
sx q[2];
rz(-0.73840529) q[2];
sx q[2];
rz(2.3389471) q[2];
rz(0.32043996) q[3];
sx q[3];
rz(-1.3405864) q[3];
sx q[3];
rz(-1.6360487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36444148) q[0];
sx q[0];
rz(-3.1138804) q[0];
sx q[0];
rz(-3.1112772) q[0];
rz(-1.3069356) q[1];
sx q[1];
rz(-2.0590643) q[1];
sx q[1];
rz(2.511715) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8960285) q[0];
sx q[0];
rz(-1.3828074) q[0];
sx q[0];
rz(2.7134368) q[0];
rz(-pi) q[1];
rz(-1.1024878) q[2];
sx q[2];
rz(-1.1716191) q[2];
sx q[2];
rz(2.1607527) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4737638) q[1];
sx q[1];
rz(-0.35982096) q[1];
sx q[1];
rz(-0.61127616) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.685018) q[3];
sx q[3];
rz(-1.1291613) q[3];
sx q[3];
rz(-0.87678443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.031875413) q[2];
sx q[2];
rz(-2.4745291) q[2];
sx q[2];
rz(-1.9943705) q[2];
rz(-2.4540497) q[3];
sx q[3];
rz(-2.5496428) q[3];
sx q[3];
rz(1.8182925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3654093) q[0];
sx q[0];
rz(-2.9054346) q[0];
sx q[0];
rz(-2.6614406) q[0];
rz(-2.7393553) q[1];
sx q[1];
rz(-1.4996585) q[1];
sx q[1];
rz(-0.38280907) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75341171) q[0];
sx q[0];
rz(-0.39330244) q[0];
sx q[0];
rz(0.9939145) q[0];
x q[1];
rz(-1.2023964) q[2];
sx q[2];
rz(-2.830626) q[2];
sx q[2];
rz(2.5045365) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.20732197) q[1];
sx q[1];
rz(-1.8658966) q[1];
sx q[1];
rz(1.1714102) q[1];
rz(1.6005101) q[3];
sx q[3];
rz(-1.374657) q[3];
sx q[3];
rz(2.8063584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.39168921) q[2];
sx q[2];
rz(-0.73035556) q[2];
sx q[2];
rz(-0.35823092) q[2];
rz(-0.41424888) q[3];
sx q[3];
rz(-1.0284938) q[3];
sx q[3];
rz(-0.39345583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8527894) q[0];
sx q[0];
rz(-1.2816592) q[0];
sx q[0];
rz(-0.417867) q[0];
rz(-1.6425543) q[1];
sx q[1];
rz(-0.42048979) q[1];
sx q[1];
rz(2.5299759) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2177416) q[0];
sx q[0];
rz(-1.799066) q[0];
sx q[0];
rz(-1.7817253) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29199227) q[2];
sx q[2];
rz(-2.9270009) q[2];
sx q[2];
rz(0.24927441) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.32595769) q[1];
sx q[1];
rz(-1.2919869) q[1];
sx q[1];
rz(1.0788979) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6526653) q[3];
sx q[3];
rz(-1.9165526) q[3];
sx q[3];
rz(0.053178259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9548107) q[2];
sx q[2];
rz(-1.6839226) q[2];
sx q[2];
rz(2.7853454) q[2];
rz(-2.6567843) q[3];
sx q[3];
rz(-1.3508513) q[3];
sx q[3];
rz(3.1118605) q[3];
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
rz(2.818882) q[0];
sx q[0];
rz(-1.1047381) q[0];
sx q[0];
rz(-1.4622965) q[0];
rz(2.7541584) q[1];
sx q[1];
rz(-2.2144364) q[1];
sx q[1];
rz(0.80518728) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96435409) q[0];
sx q[0];
rz(-2.2573009) q[0];
sx q[0];
rz(-0.80857386) q[0];
rz(-pi) q[1];
rz(-2.0207241) q[2];
sx q[2];
rz(-1.9431149) q[2];
sx q[2];
rz(-1.9228976) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3534814) q[1];
sx q[1];
rz(-2.821021) q[1];
sx q[1];
rz(-2.660257) q[1];
x q[2];
rz(-1.322454) q[3];
sx q[3];
rz(-1.9232009) q[3];
sx q[3];
rz(-0.067841522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.33599535) q[2];
sx q[2];
rz(-0.62312859) q[2];
sx q[2];
rz(1.222329) q[2];
rz(2.3313816) q[3];
sx q[3];
rz(-0.92456341) q[3];
sx q[3];
rz(-1.5789825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41351086) q[0];
sx q[0];
rz(-1.5039197) q[0];
sx q[0];
rz(1.529827) q[0];
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
rz(-2.1199216) q[3];
sx q[3];
rz(-2.3952978) q[3];
sx q[3];
rz(1.8412347) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
