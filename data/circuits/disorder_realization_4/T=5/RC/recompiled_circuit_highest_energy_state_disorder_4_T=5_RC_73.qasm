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
rz(0.29851222) q[0];
sx q[0];
rz(4.4311509) q[0];
sx q[0];
rz(9.4480954) q[0];
rz(-2.159637) q[1];
sx q[1];
rz(-2.4040931) q[1];
sx q[1];
rz(-1.4758543) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1630572) q[0];
sx q[0];
rz(-0.082162372) q[0];
sx q[0];
rz(-2.1294247) q[0];
rz(-pi) q[1];
x q[1];
rz(1.453773) q[2];
sx q[2];
rz(-2.8092779) q[2];
sx q[2];
rz(-1.8154643) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3838073) q[1];
sx q[1];
rz(-1.1769036) q[1];
sx q[1];
rz(-2.3747613) q[1];
rz(-pi) q[2];
rz(-2.2470583) q[3];
sx q[3];
rz(-1.1522376) q[3];
sx q[3];
rz(2.7194617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4618571) q[2];
sx q[2];
rz(-1.6697786) q[2];
sx q[2];
rz(2.8409345) q[2];
rz(-2.4833208) q[3];
sx q[3];
rz(-0.39977795) q[3];
sx q[3];
rz(-2.5532653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41475007) q[0];
sx q[0];
rz(-2.2756133) q[0];
sx q[0];
rz(-2.9657189) q[0];
rz(-2.580592) q[1];
sx q[1];
rz(-0.70204061) q[1];
sx q[1];
rz(-2.9452513) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6628274) q[0];
sx q[0];
rz(-1.9050373) q[0];
sx q[0];
rz(2.0408551) q[0];
rz(-pi) q[1];
rz(-1.4540599) q[2];
sx q[2];
rz(-1.4521086) q[2];
sx q[2];
rz(0.27107474) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.64164987) q[1];
sx q[1];
rz(-1.5727991) q[1];
sx q[1];
rz(2.0862122) q[1];
x q[2];
rz(-1.4370055) q[3];
sx q[3];
rz(-1.9113298) q[3];
sx q[3];
rz(1.2029369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.089513) q[2];
sx q[2];
rz(-2.5255346) q[2];
sx q[2];
rz(-1.4698131) q[2];
rz(0.52516627) q[3];
sx q[3];
rz(-1.7142121) q[3];
sx q[3];
rz(-2.4770881) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0199652) q[0];
sx q[0];
rz(-0.88931924) q[0];
sx q[0];
rz(0.63051939) q[0];
rz(-1.1446674) q[1];
sx q[1];
rz(-1.8960709) q[1];
sx q[1];
rz(-3.0847881) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5543514) q[0];
sx q[0];
rz(-1.2040753) q[0];
sx q[0];
rz(-1.1318867) q[0];
x q[1];
rz(0.47614758) q[2];
sx q[2];
rz(-2.2000929) q[2];
sx q[2];
rz(2.4193633) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2395798) q[1];
sx q[1];
rz(-0.48643349) q[1];
sx q[1];
rz(2.835266) q[1];
rz(0.2127368) q[3];
sx q[3];
rz(-1.1520583) q[3];
sx q[3];
rz(-1.1513082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.28222617) q[2];
sx q[2];
rz(-1.4713918) q[2];
sx q[2];
rz(-0.43914208) q[2];
rz(-1.7450843) q[3];
sx q[3];
rz(-1.8789995) q[3];
sx q[3];
rz(2.5631574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0466995) q[0];
sx q[0];
rz(-0.72143227) q[0];
sx q[0];
rz(-1.0149581) q[0];
rz(0.062648423) q[1];
sx q[1];
rz(-0.95476127) q[1];
sx q[1];
rz(2.8573724) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6054404) q[0];
sx q[0];
rz(-1.221719) q[0];
sx q[0];
rz(2.6830868) q[0];
rz(-1.4309056) q[2];
sx q[2];
rz(-1.0770742) q[2];
sx q[2];
rz(1.7715275) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.758421) q[1];
sx q[1];
rz(-1.212442) q[1];
sx q[1];
rz(1.5835632) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8241513) q[3];
sx q[3];
rz(-1.5572963) q[3];
sx q[3];
rz(2.6576192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8404428) q[2];
sx q[2];
rz(-2.132685) q[2];
sx q[2];
rz(-0.76048771) q[2];
rz(2.8216951) q[3];
sx q[3];
rz(-2.0127998) q[3];
sx q[3];
rz(2.1130051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8496721) q[0];
sx q[0];
rz(-2.5909162) q[0];
sx q[0];
rz(-0.40529761) q[0];
rz(-2.6536476) q[1];
sx q[1];
rz(-1.1064203) q[1];
sx q[1];
rz(0.36283666) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0340509) q[0];
sx q[0];
rz(-1.1404317) q[0];
sx q[0];
rz(-1.2314046) q[0];
x q[1];
rz(-2.5561781) q[2];
sx q[2];
rz(-1.2512794) q[2];
sx q[2];
rz(-1.320517) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2606973) q[1];
sx q[1];
rz(-2.915307) q[1];
sx q[1];
rz(0.21842167) q[1];
rz(-0.69546206) q[3];
sx q[3];
rz(-1.3331279) q[3];
sx q[3];
rz(-0.15791751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.001658) q[2];
sx q[2];
rz(-1.1589061) q[2];
sx q[2];
rz(0.72251594) q[2];
rz(-2.6288988) q[3];
sx q[3];
rz(-2.2721458) q[3];
sx q[3];
rz(1.9374013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8743415) q[0];
sx q[0];
rz(-2.5047472) q[0];
sx q[0];
rz(-0.27477086) q[0];
rz(-0.35708669) q[1];
sx q[1];
rz(-1.6386702) q[1];
sx q[1];
rz(-1.1579317) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0561075) q[0];
sx q[0];
rz(-2.4153313) q[0];
sx q[0];
rz(1.1994491) q[0];
x q[1];
rz(0.59717565) q[2];
sx q[2];
rz(-0.9998601) q[2];
sx q[2];
rz(-1.2025637) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9951862) q[1];
sx q[1];
rz(-0.9764834) q[1];
sx q[1];
rz(1.8176778) q[1];
x q[2];
rz(2.255106) q[3];
sx q[3];
rz(-2.2456256) q[3];
sx q[3];
rz(1.0533028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.57924119) q[2];
sx q[2];
rz(-1.7513821) q[2];
sx q[2];
rz(-0.58970279) q[2];
rz(1.3127182) q[3];
sx q[3];
rz(-1.4785942) q[3];
sx q[3];
rz(-2.5780799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24484672) q[0];
sx q[0];
rz(-0.99269301) q[0];
sx q[0];
rz(-0.27031159) q[0];
rz(-2.7145794) q[1];
sx q[1];
rz(-1.3753563) q[1];
sx q[1];
rz(-2.7332773) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9664551) q[0];
sx q[0];
rz(-1.9388559) q[0];
sx q[0];
rz(-1.4428663) q[0];
x q[1];
rz(-1.4502268) q[2];
sx q[2];
rz(-0.70115108) q[2];
sx q[2];
rz(2.1599959) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9198614) q[1];
sx q[1];
rz(-1.9834922) q[1];
sx q[1];
rz(-2.2151715) q[1];
rz(1.6736843) q[3];
sx q[3];
rz(-0.85082084) q[3];
sx q[3];
rz(-2.306769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2030877) q[2];
sx q[2];
rz(-2.0577343) q[2];
sx q[2];
rz(-0.92424029) q[2];
rz(-2.3067394) q[3];
sx q[3];
rz(-0.82930851) q[3];
sx q[3];
rz(0.99669641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95558178) q[0];
sx q[0];
rz(-0.38014933) q[0];
sx q[0];
rz(0.36744776) q[0];
rz(0.5683178) q[1];
sx q[1];
rz(-1.808993) q[1];
sx q[1];
rz(0.41499358) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6870553) q[0];
sx q[0];
rz(-0.53249411) q[0];
sx q[0];
rz(-2.103602) q[0];
x q[1];
rz(2.6122323) q[2];
sx q[2];
rz(-1.0870013) q[2];
sx q[2];
rz(2.8279378) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7447978) q[1];
sx q[1];
rz(-0.84229453) q[1];
sx q[1];
rz(-0.70835967) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0051963768) q[3];
sx q[3];
rz(-1.2250161) q[3];
sx q[3];
rz(-0.41692641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.35855287) q[2];
sx q[2];
rz(-0.8873322) q[2];
sx q[2];
rz(2.3084194) q[2];
rz(0.35946515) q[3];
sx q[3];
rz(-1.0901674) q[3];
sx q[3];
rz(-1.6370714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1737162) q[0];
sx q[0];
rz(-2.739527) q[0];
sx q[0];
rz(3.1267082) q[0];
rz(-1.124565) q[1];
sx q[1];
rz(-2.896307) q[1];
sx q[1];
rz(-3.0344149) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5456148) q[0];
sx q[0];
rz(-2.2920485) q[0];
sx q[0];
rz(0.19907339) q[0];
rz(0.46469073) q[2];
sx q[2];
rz(-0.87848262) q[2];
sx q[2];
rz(-1.9141045) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8228777) q[1];
sx q[1];
rz(-1.7194304) q[1];
sx q[1];
rz(0.48245247) q[1];
x q[2];
rz(1.393287) q[3];
sx q[3];
rz(-1.313398) q[3];
sx q[3];
rz(2.7696115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5643481) q[2];
sx q[2];
rz(-1.8337092) q[2];
sx q[2];
rz(0.85961071) q[2];
rz(0.25935069) q[3];
sx q[3];
rz(-1.3602463) q[3];
sx q[3];
rz(-2.7249961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6421826) q[0];
sx q[0];
rz(-3.0152617) q[0];
sx q[0];
rz(-1.1580178) q[0];
rz(-2.6462789) q[1];
sx q[1];
rz(-1.9751578) q[1];
sx q[1];
rz(-0.29702979) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7796233) q[0];
sx q[0];
rz(-1.8898148) q[0];
sx q[0];
rz(-3.1373128) q[0];
rz(-pi) q[1];
rz(-3.0115602) q[2];
sx q[2];
rz(-0.30270019) q[2];
sx q[2];
rz(-1.5695733) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.37849879) q[1];
sx q[1];
rz(-1.2796254) q[1];
sx q[1];
rz(-1.4062455) q[1];
x q[2];
rz(-0.90822753) q[3];
sx q[3];
rz(-1.5836827) q[3];
sx q[3];
rz(2.9290309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.872252) q[2];
sx q[2];
rz(-0.38186914) q[2];
sx q[2];
rz(-2.2807109) q[2];
rz(0.46368972) q[3];
sx q[3];
rz(-0.90353614) q[3];
sx q[3];
rz(-1.1369107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9853482) q[0];
sx q[0];
rz(-1.5679659) q[0];
sx q[0];
rz(1.1563942) q[0];
rz(-1.8078177) q[1];
sx q[1];
rz(-1.6009686) q[1];
sx q[1];
rz(0.70645465) q[1];
rz(2.4598224) q[2];
sx q[2];
rz(-0.78777704) q[2];
sx q[2];
rz(-2.4252979) q[2];
rz(-0.49336962) q[3];
sx q[3];
rz(-1.9153638) q[3];
sx q[3];
rz(2.3079688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
