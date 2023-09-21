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
rz(-2.535948) q[1];
sx q[1];
rz(-0.65094596) q[1];
sx q[1];
rz(-2.5075066) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2274322) q[0];
sx q[0];
rz(-2.3254546) q[0];
sx q[0];
rz(-0.039365191) q[0];
x q[1];
rz(-3.0787266) q[2];
sx q[2];
rz(-0.91744423) q[2];
sx q[2];
rz(-1.1970929) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.69971985) q[1];
sx q[1];
rz(-2.3708214) q[1];
sx q[1];
rz(2.2295879) q[1];
rz(-pi) q[2];
rz(-0.4001873) q[3];
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
rz(1.8566711) q[3];
sx q[3];
rz(-1.4611171) q[3];
sx q[3];
rz(2.8485956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15853515) q[0];
sx q[0];
rz(-1.8479269) q[0];
sx q[0];
rz(0.43757004) q[0];
rz(-2.5105387) q[1];
sx q[1];
rz(-2.8129306) q[1];
sx q[1];
rz(-0.22110573) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32199931) q[0];
sx q[0];
rz(-2.0237192) q[0];
sx q[0];
rz(0.20690147) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0304893) q[2];
sx q[2];
rz(-1.9176033) q[2];
sx q[2];
rz(-0.32804104) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7213388) q[1];
sx q[1];
rz(-1.4495279) q[1];
sx q[1];
rz(-2.6810357) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9651871) q[3];
sx q[3];
rz(-0.50900148) q[3];
sx q[3];
rz(0.41364663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.21800403) q[2];
sx q[2];
rz(-1.4322586) q[2];
sx q[2];
rz(-2.5533) q[2];
rz(-0.44899392) q[3];
sx q[3];
rz(-0.42789999) q[3];
sx q[3];
rz(2.0935521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5730729) q[0];
sx q[0];
rz(-0.097671106) q[0];
sx q[0];
rz(-2.1411238) q[0];
rz(-2.4160642) q[1];
sx q[1];
rz(-2.0670481) q[1];
sx q[1];
rz(0.75769889) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0453148) q[0];
sx q[0];
rz(-2.0925539) q[0];
sx q[0];
rz(1.9406712) q[0];
rz(0.65501113) q[2];
sx q[2];
rz(-2.8821324) q[2];
sx q[2];
rz(0.26817817) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.12049) q[1];
sx q[1];
rz(-2.0334525) q[1];
sx q[1];
rz(-0.91336577) q[1];
x q[2];
rz(0.63852255) q[3];
sx q[3];
rz(-1.4811852) q[3];
sx q[3];
rz(1.5479969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9329325) q[2];
sx q[2];
rz(-0.22298403) q[2];
sx q[2];
rz(0.83646742) q[2];
rz(1.4423192) q[3];
sx q[3];
rz(-1.6675555) q[3];
sx q[3];
rz(-2.9377655) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02012415) q[0];
sx q[0];
rz(-0.080105372) q[0];
sx q[0];
rz(2.6304723) q[0];
rz(0.077443667) q[1];
sx q[1];
rz(-0.42235342) q[1];
sx q[1];
rz(1.3285332) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3805566) q[0];
sx q[0];
rz(-1.3576926) q[0];
sx q[0];
rz(-2.9831432) q[0];
x q[1];
rz(-2.0650234) q[2];
sx q[2];
rz(-0.87072125) q[2];
sx q[2];
rz(-1.9698524) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7698344) q[1];
sx q[1];
rz(-1.8375988) q[1];
sx q[1];
rz(3.1410602) q[1];
x q[2];
rz(-2.02416) q[3];
sx q[3];
rz(-0.24586596) q[3];
sx q[3];
rz(0.21594957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0478583) q[2];
sx q[2];
rz(-1.9231984) q[2];
sx q[2];
rz(1.1070586) q[2];
rz(0.47248653) q[3];
sx q[3];
rz(-1.5269591) q[3];
sx q[3];
rz(-2.2842177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040314019) q[0];
sx q[0];
rz(-2.2450228) q[0];
sx q[0];
rz(2.8033946) q[0];
rz(-1.2942554) q[1];
sx q[1];
rz(-1.5599524) q[1];
sx q[1];
rz(-0.50450528) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5915247) q[0];
sx q[0];
rz(-1.803777) q[0];
sx q[0];
rz(-2.1165119) q[0];
rz(-pi) q[1];
rz(0.46576969) q[2];
sx q[2];
rz(-1.2301187) q[2];
sx q[2];
rz(-2.6054232) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.075715847) q[1];
sx q[1];
rz(-1.8738235) q[1];
sx q[1];
rz(2.8916343) q[1];
x q[2];
rz(2.3584064) q[3];
sx q[3];
rz(-0.83055701) q[3];
sx q[3];
rz(-0.98379788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.111104) q[2];
sx q[2];
rz(-0.8384853) q[2];
sx q[2];
rz(1.4939235) q[2];
rz(1.6882287) q[3];
sx q[3];
rz(-0.93795347) q[3];
sx q[3];
rz(1.7593613) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9313653) q[0];
sx q[0];
rz(-2.2670822) q[0];
sx q[0];
rz(-1.0908303) q[0];
rz(-2.6018654) q[1];
sx q[1];
rz(-1.9878186) q[1];
sx q[1];
rz(2.9528023) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.328619) q[0];
sx q[0];
rz(-0.44437528) q[0];
sx q[0];
rz(0.61291738) q[0];
x q[1];
rz(0.10411711) q[2];
sx q[2];
rz(-1.3008442) q[2];
sx q[2];
rz(-1.5013258) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1660359) q[1];
sx q[1];
rz(-1.3284725) q[1];
sx q[1];
rz(-1.5029328) q[1];
x q[2];
rz(-2.4981899) q[3];
sx q[3];
rz(-0.2989558) q[3];
sx q[3];
rz(1.1645731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9699036) q[2];
sx q[2];
rz(-2.1128555) q[2];
sx q[2];
rz(1.3151273) q[2];
rz(1.6361489) q[3];
sx q[3];
rz(-2.7841778) q[3];
sx q[3];
rz(-1.858254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0514907) q[0];
sx q[0];
rz(-2.3110456) q[0];
sx q[0];
rz(0.32399696) q[0];
rz(-1.8404768) q[1];
sx q[1];
rz(-0.83530656) q[1];
sx q[1];
rz(1.3791929) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3430816) q[0];
sx q[0];
rz(-1.8146975) q[0];
sx q[0];
rz(0.084246158) q[0];
rz(-pi) q[1];
rz(-2.4263779) q[2];
sx q[2];
rz(-2.5589802) q[2];
sx q[2];
rz(1.6974534) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0849689) q[1];
sx q[1];
rz(-1.2722204) q[1];
sx q[1];
rz(1.741239) q[1];
rz(2.2961388) q[3];
sx q[3];
rz(-0.77973706) q[3];
sx q[3];
rz(-0.23507915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.032701187) q[2];
sx q[2];
rz(-1.4543433) q[2];
sx q[2];
rz(-2.1984055) q[2];
rz(2.8105248) q[3];
sx q[3];
rz(-1.3676684) q[3];
sx q[3];
rz(0.29512063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11809764) q[0];
sx q[0];
rz(-2.5546615) q[0];
sx q[0];
rz(2.3773637) q[0];
rz(-3.0006192) q[1];
sx q[1];
rz(-0.39607513) q[1];
sx q[1];
rz(-1.8483298) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9613567) q[0];
sx q[0];
rz(-1.1399674) q[0];
sx q[0];
rz(-0.21813099) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60111945) q[2];
sx q[2];
rz(-1.6496611) q[2];
sx q[2];
rz(-2.2476946) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.075899374) q[1];
sx q[1];
rz(-1.7095487) q[1];
sx q[1];
rz(-1.1636415) q[1];
rz(-pi) q[2];
rz(-0.44869081) q[3];
sx q[3];
rz(-0.43827) q[3];
sx q[3];
rz(1.1718307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3747037) q[2];
sx q[2];
rz(-2.1024487) q[2];
sx q[2];
rz(0.58132201) q[2];
rz(2.2733722) q[3];
sx q[3];
rz(-0.41728443) q[3];
sx q[3];
rz(1.8301615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54810369) q[0];
sx q[0];
rz(-0.5797525) q[0];
sx q[0];
rz(-1.2506437) q[0];
rz(2.1447694) q[1];
sx q[1];
rz(-0.88368982) q[1];
sx q[1];
rz(-1.2876127) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2851965) q[0];
sx q[0];
rz(-2.6515793) q[0];
sx q[0];
rz(-1.6243837) q[0];
x q[1];
rz(-2.8224045) q[2];
sx q[2];
rz(-0.50779283) q[2];
sx q[2];
rz(-2.9920981) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7848283) q[1];
sx q[1];
rz(-0.67912662) q[1];
sx q[1];
rz(-1.9041512) q[1];
rz(-0.82448126) q[3];
sx q[3];
rz(-1.5319676) q[3];
sx q[3];
rz(-0.83394921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.320497) q[2];
sx q[2];
rz(-0.38791502) q[2];
sx q[2];
rz(1.256475) q[2];
rz(0.16658941) q[3];
sx q[3];
rz(-1.5766671) q[3];
sx q[3];
rz(2.0690209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2607516) q[0];
sx q[0];
rz(-2.3374225) q[0];
sx q[0];
rz(-2.8163731) q[0];
rz(-1.1351769) q[1];
sx q[1];
rz(-2.4644641) q[1];
sx q[1];
rz(-0.36718711) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4627535) q[0];
sx q[0];
rz(-0.05719962) q[0];
sx q[0];
rz(-2.5662533) q[0];
rz(-pi) q[1];
rz(0.37819241) q[2];
sx q[2];
rz(-0.66236712) q[2];
sx q[2];
rz(-2.3364002) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2314184) q[1];
sx q[1];
rz(-1.3712198) q[1];
sx q[1];
rz(2.477596) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1807947) q[3];
sx q[3];
rz(-0.73878091) q[3];
sx q[3];
rz(1.9566386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8964768) q[2];
sx q[2];
rz(-2.3325236) q[2];
sx q[2];
rz(2.0754576) q[2];
rz(-0.079244763) q[3];
sx q[3];
rz(-2.5388122) q[3];
sx q[3];
rz(1.4762896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8463678) q[0];
sx q[0];
rz(-1.5867148) q[0];
sx q[0];
rz(-1.5750194) q[0];
rz(-2.8425343) q[1];
sx q[1];
rz(-2.538264) q[1];
sx q[1];
rz(-2.6187142) q[1];
rz(1.4916186) q[2];
sx q[2];
rz(-0.82293246) q[2];
sx q[2];
rz(-0.056597829) q[2];
rz(-2.9885837) q[3];
sx q[3];
rz(-0.85481337) q[3];
sx q[3];
rz(-0.92683642) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
