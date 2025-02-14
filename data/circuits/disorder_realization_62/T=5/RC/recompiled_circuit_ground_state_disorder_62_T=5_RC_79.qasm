OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5247076) q[0];
sx q[0];
rz(-2.4029713) q[0];
sx q[0];
rz(3.024616) q[0];
rz(0.95247954) q[1];
sx q[1];
rz(-1.6713961) q[1];
sx q[1];
rz(-2.267946) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0067209816) q[0];
sx q[0];
rz(-0.76113478) q[0];
sx q[0];
rz(-1.8675682) q[0];
rz(-pi) q[1];
rz(-1.123465) q[2];
sx q[2];
rz(-1.8904933) q[2];
sx q[2];
rz(1.5716219) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4651771) q[1];
sx q[1];
rz(-2.5897674) q[1];
sx q[1];
rz(-1.4613749) q[1];
x q[2];
rz(1.1659032) q[3];
sx q[3];
rz(-1.3171853) q[3];
sx q[3];
rz(3.1133127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9690507) q[2];
sx q[2];
rz(-1.8092864) q[2];
sx q[2];
rz(1.8633899) q[2];
rz(1.5365907) q[3];
sx q[3];
rz(-1.2383818) q[3];
sx q[3];
rz(2.8240375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69739598) q[0];
sx q[0];
rz(-0.2187271) q[0];
sx q[0];
rz(0.87795192) q[0];
rz(0.76389337) q[1];
sx q[1];
rz(-2.0593675) q[1];
sx q[1];
rz(2.0108932) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5661086) q[0];
sx q[0];
rz(-2.974925) q[0];
sx q[0];
rz(-1.9420293) q[0];
rz(-1.892852) q[2];
sx q[2];
rz(-0.45834228) q[2];
sx q[2];
rz(-0.94142585) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.69399) q[1];
sx q[1];
rz(-2.388235) q[1];
sx q[1];
rz(-1.5521469) q[1];
rz(-2.9933287) q[3];
sx q[3];
rz(-0.63126031) q[3];
sx q[3];
rz(0.12693044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6905602) q[2];
sx q[2];
rz(-1.6310383) q[2];
sx q[2];
rz(-0.65563273) q[2];
rz(-1.2304652) q[3];
sx q[3];
rz(-2.1793607) q[3];
sx q[3];
rz(-0.83278304) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3967628) q[0];
sx q[0];
rz(-1.1711045) q[0];
sx q[0];
rz(0.50286621) q[0];
rz(-0.14547959) q[1];
sx q[1];
rz(-1.611404) q[1];
sx q[1];
rz(-1.9810289) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6737719) q[0];
sx q[0];
rz(-1.880411) q[0];
sx q[0];
rz(-1.2542115) q[0];
rz(1.4394748) q[2];
sx q[2];
rz(-1.223837) q[2];
sx q[2];
rz(1.0441213) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.33659157) q[1];
sx q[1];
rz(-2.1865926) q[1];
sx q[1];
rz(-2.6686882) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3011906) q[3];
sx q[3];
rz(-1.9072215) q[3];
sx q[3];
rz(-2.2158282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1955371) q[2];
sx q[2];
rz(-1.035752) q[2];
sx q[2];
rz(-2.3717086) q[2];
rz(-0.59512538) q[3];
sx q[3];
rz(-0.49134058) q[3];
sx q[3];
rz(-2.6814521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50764099) q[0];
sx q[0];
rz(-2.0946298) q[0];
sx q[0];
rz(1.7858343) q[0];
rz(-2.5178364) q[1];
sx q[1];
rz(-1.535894) q[1];
sx q[1];
rz(0.49130586) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3771916) q[0];
sx q[0];
rz(-2.1024414) q[0];
sx q[0];
rz(1.3926943) q[0];
x q[1];
rz(-2.5042438) q[2];
sx q[2];
rz(-2.7985811) q[2];
sx q[2];
rz(-2.1528139) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2249105) q[1];
sx q[1];
rz(-2.2631117) q[1];
sx q[1];
rz(0.40596753) q[1];
x q[2];
rz(1.0503429) q[3];
sx q[3];
rz(-2.1698275) q[3];
sx q[3];
rz(1.1831621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9451311) q[2];
sx q[2];
rz(-1.6170231) q[2];
sx q[2];
rz(-0.54836908) q[2];
rz(-1.3301001) q[3];
sx q[3];
rz(-2.8434704) q[3];
sx q[3];
rz(0.33838457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86692989) q[0];
sx q[0];
rz(-0.45329705) q[0];
sx q[0];
rz(1.0484265) q[0];
rz(0.10920814) q[1];
sx q[1];
rz(-1.5501153) q[1];
sx q[1];
rz(1.106326) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7503742) q[0];
sx q[0];
rz(-2.8958909) q[0];
sx q[0];
rz(2.0294163) q[0];
rz(1.1961522) q[2];
sx q[2];
rz(-0.76065791) q[2];
sx q[2];
rz(0.51860318) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2611003) q[1];
sx q[1];
rz(-1.2478117) q[1];
sx q[1];
rz(2.9862981) q[1];
rz(-2.2436687) q[3];
sx q[3];
rz(-0.37386383) q[3];
sx q[3];
rz(-1.7927235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.76935524) q[2];
sx q[2];
rz(-1.8753588) q[2];
sx q[2];
rz(-0.20277578) q[2];
rz(0.77962223) q[3];
sx q[3];
rz(-2.0519665) q[3];
sx q[3];
rz(0.37105086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5146273) q[0];
sx q[0];
rz(-3.096464) q[0];
sx q[0];
rz(-0.87920642) q[0];
rz(-2.5458287) q[1];
sx q[1];
rz(-0.87409449) q[1];
sx q[1];
rz(-1.2557868) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6539772) q[0];
sx q[0];
rz(-0.10627986) q[0];
sx q[0];
rz(2.6208861) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57298341) q[2];
sx q[2];
rz(-1.8610916) q[2];
sx q[2];
rz(-0.03183768) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.75460282) q[1];
sx q[1];
rz(-2.5856254) q[1];
sx q[1];
rz(2.703859) q[1];
x q[2];
rz(-0.34481315) q[3];
sx q[3];
rz(-0.64467829) q[3];
sx q[3];
rz(-0.0078474069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6229728) q[2];
sx q[2];
rz(-2.8039248) q[2];
sx q[2];
rz(-1.6597623) q[2];
rz(1.6984113) q[3];
sx q[3];
rz(-1.7549763) q[3];
sx q[3];
rz(1.1186918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54565322) q[0];
sx q[0];
rz(-2.312199) q[0];
sx q[0];
rz(-0.19350061) q[0];
rz(0.045348383) q[1];
sx q[1];
rz(-1.9769316) q[1];
sx q[1];
rz(0.68101105) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15475965) q[0];
sx q[0];
rz(-1.7333507) q[0];
sx q[0];
rz(-1.6587371) q[0];
rz(-1.1521336) q[2];
sx q[2];
rz(-2.5730926) q[2];
sx q[2];
rz(-0.9886888) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4894708) q[1];
sx q[1];
rz(-0.014111405) q[1];
sx q[1];
rz(1.1443787) q[1];
rz(0.29905921) q[3];
sx q[3];
rz(-0.402723) q[3];
sx q[3];
rz(-2.2880768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.27651522) q[2];
sx q[2];
rz(-2.4962208) q[2];
sx q[2];
rz(1.7934249) q[2];
rz(1.3639785) q[3];
sx q[3];
rz(-1.8230702) q[3];
sx q[3];
rz(1.9510423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6337223) q[0];
sx q[0];
rz(-2.3220334) q[0];
sx q[0];
rz(0.67935294) q[0];
rz(-2.9044115) q[1];
sx q[1];
rz(-0.91865426) q[1];
sx q[1];
rz(-1.0666581) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53700457) q[0];
sx q[0];
rz(-2.6312772) q[0];
sx q[0];
rz(-3.133581) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1751498) q[2];
sx q[2];
rz(-1.4928515) q[2];
sx q[2];
rz(-1.9029531) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.28344) q[1];
sx q[1];
rz(-1.6951218) q[1];
sx q[1];
rz(2.0945626) q[1];
rz(1.1447843) q[3];
sx q[3];
rz(-1.2698114) q[3];
sx q[3];
rz(2.3097484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8827768) q[2];
sx q[2];
rz(-1.2315742) q[2];
sx q[2];
rz(0.41137496) q[2];
rz(1.2808293) q[3];
sx q[3];
rz(-2.4954093) q[3];
sx q[3];
rz(-1.0584186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48788747) q[0];
sx q[0];
rz(-1.7937086) q[0];
sx q[0];
rz(-2.2109798) q[0];
rz(0.46512428) q[1];
sx q[1];
rz(-2.5587176) q[1];
sx q[1];
rz(1.7238269) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0865308) q[0];
sx q[0];
rz(-1.3320005) q[0];
sx q[0];
rz(-0.60542528) q[0];
rz(0.82724656) q[2];
sx q[2];
rz(-2.3066024) q[2];
sx q[2];
rz(-2.5967251) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5807996) q[1];
sx q[1];
rz(-1.8617814) q[1];
sx q[1];
rz(2.7853904) q[1];
rz(-3.0259589) q[3];
sx q[3];
rz(-0.99471417) q[3];
sx q[3];
rz(-0.86365684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.96665367) q[2];
sx q[2];
rz(-0.83028364) q[2];
sx q[2];
rz(-1.9261599) q[2];
rz(0.71839607) q[3];
sx q[3];
rz(-1.6738439) q[3];
sx q[3];
rz(2.4992656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65828085) q[0];
sx q[0];
rz(-2.719306) q[0];
sx q[0];
rz(1.8102113) q[0];
rz(-0.70679682) q[1];
sx q[1];
rz(-1.4247318) q[1];
sx q[1];
rz(1.4250863) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0646149) q[0];
sx q[0];
rz(-2.0138903) q[0];
sx q[0];
rz(-0.71378543) q[0];
x q[1];
rz(-3.0767137) q[2];
sx q[2];
rz(-2.2264495) q[2];
sx q[2];
rz(-0.82198373) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6687732) q[1];
sx q[1];
rz(-2.4309078) q[1];
sx q[1];
rz(-2.4345934) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44300191) q[3];
sx q[3];
rz(-0.85369977) q[3];
sx q[3];
rz(2.9505961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4786272) q[2];
sx q[2];
rz(-2.4586283) q[2];
sx q[2];
rz(1.9662201) q[2];
rz(3.1183682) q[3];
sx q[3];
rz(-0.5718137) q[3];
sx q[3];
rz(2.8779252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57446734) q[0];
sx q[0];
rz(-2.1463285) q[0];
sx q[0];
rz(1.2405201) q[0];
rz(2.4868838) q[1];
sx q[1];
rz(-1.2798825) q[1];
sx q[1];
rz(1.9725694) q[1];
rz(2.2562182) q[2];
sx q[2];
rz(-0.17687951) q[2];
sx q[2];
rz(2.542873) q[2];
rz(-2.2795879) q[3];
sx q[3];
rz(-1.0222407) q[3];
sx q[3];
rz(0.23289451) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
