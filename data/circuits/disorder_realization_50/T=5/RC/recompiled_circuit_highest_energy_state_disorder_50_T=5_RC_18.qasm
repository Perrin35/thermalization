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
rz(0.69831508) q[0];
sx q[0];
rz(-0.3787711) q[0];
sx q[0];
rz(-1.9842499) q[0];
rz(10.209822) q[1];
sx q[1];
rz(0.68612376) q[1];
sx q[1];
rz(2.6148028) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0924659) q[0];
sx q[0];
rz(-1.2727588) q[0];
sx q[0];
rz(-2.1597693) q[0];
rz(0.77157048) q[2];
sx q[2];
rz(-2.2419746) q[2];
sx q[2];
rz(-3.0085466) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96054582) q[1];
sx q[1];
rz(-1.6072453) q[1];
sx q[1];
rz(-0.44238449) q[1];
rz(-pi) q[2];
rz(1.6185226) q[3];
sx q[3];
rz(-2.5807947) q[3];
sx q[3];
rz(0.58855295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.779125) q[2];
sx q[2];
rz(-2.1802826) q[2];
sx q[2];
rz(0.15164068) q[2];
rz(-0.14196299) q[3];
sx q[3];
rz(-1.4700593) q[3];
sx q[3];
rz(-1.1375455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66459429) q[0];
sx q[0];
rz(-2.4095896) q[0];
sx q[0];
rz(2.3240996) q[0];
rz(-3.0290161) q[1];
sx q[1];
rz(-1.6855626) q[1];
sx q[1];
rz(-2.2419825) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7821723) q[0];
sx q[0];
rz(-0.36672584) q[0];
sx q[0];
rz(0.72215778) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1915273) q[2];
sx q[2];
rz(-2.143444) q[2];
sx q[2];
rz(0.060163035) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6625329) q[1];
sx q[1];
rz(-2.3182437) q[1];
sx q[1];
rz(-1.6475299) q[1];
rz(-pi) q[2];
rz(-2.6255304) q[3];
sx q[3];
rz(-0.96092292) q[3];
sx q[3];
rz(-0.42760951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4704935) q[2];
sx q[2];
rz(-1.495139) q[2];
sx q[2];
rz(-2.1136843) q[2];
rz(0.73728621) q[3];
sx q[3];
rz(-1.6643915) q[3];
sx q[3];
rz(3.1173053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0055493) q[0];
sx q[0];
rz(-2.5856954) q[0];
sx q[0];
rz(1.9480202) q[0];
rz(-1.8434803) q[1];
sx q[1];
rz(-1.6622512) q[1];
sx q[1];
rz(-1.6568291) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1910716) q[0];
sx q[0];
rz(-1.6170039) q[0];
sx q[0];
rz(-0.025006983) q[0];
rz(3.0151691) q[2];
sx q[2];
rz(-1.3609386) q[2];
sx q[2];
rz(-1.017638) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.79968184) q[1];
sx q[1];
rz(-1.4355625) q[1];
sx q[1];
rz(2.7970201) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0590634) q[3];
sx q[3];
rz(-1.2090873) q[3];
sx q[3];
rz(-1.495273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9803311) q[2];
sx q[2];
rz(-1.9183466) q[2];
sx q[2];
rz(0.70183357) q[2];
rz(-0.069132239) q[3];
sx q[3];
rz(-2.2528503) q[3];
sx q[3];
rz(-2.4338636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.055701) q[0];
sx q[0];
rz(-1.1393071) q[0];
sx q[0];
rz(2.5162146) q[0];
rz(-1.0944132) q[1];
sx q[1];
rz(-1.6182599) q[1];
sx q[1];
rz(1.3166924) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65526456) q[0];
sx q[0];
rz(-1.0287971) q[0];
sx q[0];
rz(-0.16714759) q[0];
rz(-pi) q[1];
rz(-2.8595692) q[2];
sx q[2];
rz(-0.33581734) q[2];
sx q[2];
rz(-1.6086325) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.81165555) q[1];
sx q[1];
rz(-1.3906758) q[1];
sx q[1];
rz(0.7094769) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1573477) q[3];
sx q[3];
rz(-0.86185019) q[3];
sx q[3];
rz(-1.1701442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7350498) q[2];
sx q[2];
rz(-2.4989765) q[2];
sx q[2];
rz(-3.1114846) q[2];
rz(-0.41664577) q[3];
sx q[3];
rz(-1.7130339) q[3];
sx q[3];
rz(3.0863975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77867126) q[0];
sx q[0];
rz(-1.8632977) q[0];
sx q[0];
rz(1.2177421) q[0];
rz(0.22449224) q[1];
sx q[1];
rz(-1.6480564) q[1];
sx q[1];
rz(-2.7405558) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5767956) q[0];
sx q[0];
rz(-2.5677997) q[0];
sx q[0];
rz(-2.3032715) q[0];
x q[1];
rz(-0.46868344) q[2];
sx q[2];
rz(-2.4309553) q[2];
sx q[2];
rz(0.24220322) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4523938) q[1];
sx q[1];
rz(-1.3340063) q[1];
sx q[1];
rz(-1.6195606) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55373904) q[3];
sx q[3];
rz(-2.9512292) q[3];
sx q[3];
rz(-2.7680264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2991221) q[2];
sx q[2];
rz(-1.9186019) q[2];
sx q[2];
rz(2.6341338) q[2];
rz(-0.92042813) q[3];
sx q[3];
rz(-2.7046552) q[3];
sx q[3];
rz(-2.9612655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64666635) q[0];
sx q[0];
rz(-0.05519069) q[0];
sx q[0];
rz(2.1221509) q[0];
rz(1.9693718) q[1];
sx q[1];
rz(-1.1727138) q[1];
sx q[1];
rz(1.2219465) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2766314) q[0];
sx q[0];
rz(-0.91827938) q[0];
sx q[0];
rz(-0.44256532) q[0];
rz(-pi) q[1];
rz(-1.5084294) q[2];
sx q[2];
rz(-2.6259239) q[2];
sx q[2];
rz(-0.0002278681) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3675569) q[1];
sx q[1];
rz(-1.8453215) q[1];
sx q[1];
rz(-0.56660272) q[1];
x q[2];
rz(-2.1453247) q[3];
sx q[3];
rz(-1.822851) q[3];
sx q[3];
rz(-1.1561405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4814066) q[2];
sx q[2];
rz(-0.6654827) q[2];
sx q[2];
rz(1.0931724) q[2];
rz(-0.53705755) q[3];
sx q[3];
rz(-1.4635181) q[3];
sx q[3];
rz(2.4721036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2375803) q[0];
sx q[0];
rz(-2.8222988) q[0];
sx q[0];
rz(-1.681666) q[0];
rz(1.6319252) q[1];
sx q[1];
rz(-2.5131112) q[1];
sx q[1];
rz(0.07930886) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96937856) q[0];
sx q[0];
rz(-1.8076118) q[0];
sx q[0];
rz(-3.0585312) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4334636) q[2];
sx q[2];
rz(-1.3690565) q[2];
sx q[2];
rz(1.2936426) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7056337) q[1];
sx q[1];
rz(-0.50928947) q[1];
sx q[1];
rz(2.5540256) q[1];
rz(-pi) q[2];
rz(2.9856624) q[3];
sx q[3];
rz(-0.72375384) q[3];
sx q[3];
rz(-0.93800217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1377533) q[2];
sx q[2];
rz(-1.1649818) q[2];
sx q[2];
rz(-2.7286781) q[2];
rz(-1.2952992) q[3];
sx q[3];
rz(-0.92419878) q[3];
sx q[3];
rz(-2.7220461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.15692784) q[0];
sx q[0];
rz(-2.4913737) q[0];
sx q[0];
rz(0.28537634) q[0];
rz(-3.0889619) q[1];
sx q[1];
rz(-1.4080518) q[1];
sx q[1];
rz(0.1836798) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3899468) q[0];
sx q[0];
rz(-3.0245298) q[0];
sx q[0];
rz(-1.2617017) q[0];
rz(-pi) q[1];
rz(-1.7247612) q[2];
sx q[2];
rz(-1.2885258) q[2];
sx q[2];
rz(-0.33253129) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0252776) q[1];
sx q[1];
rz(-2.5417788) q[1];
sx q[1];
rz(0.13529899) q[1];
x q[2];
rz(-2.7643572) q[3];
sx q[3];
rz(-2.4580015) q[3];
sx q[3];
rz(0.93631668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8949184) q[2];
sx q[2];
rz(-1.7343212) q[2];
sx q[2];
rz(1.0651917) q[2];
rz(1.6861606) q[3];
sx q[3];
rz(-1.287241) q[3];
sx q[3];
rz(-2.5559032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10373779) q[0];
sx q[0];
rz(-2.3268564) q[0];
sx q[0];
rz(-1.5392342) q[0];
rz(0.94929758) q[1];
sx q[1];
rz(-2.0930591) q[1];
sx q[1];
rz(-1.4303713) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8447782) q[0];
sx q[0];
rz(-2.0955293) q[0];
sx q[0];
rz(1.0206166) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4865506) q[2];
sx q[2];
rz(-0.73053321) q[2];
sx q[2];
rz(-0.73405311) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1518703) q[1];
sx q[1];
rz(-0.29726004) q[1];
sx q[1];
rz(-0.10135915) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28258459) q[3];
sx q[3];
rz(-1.1573829) q[3];
sx q[3];
rz(1.2334718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0091693) q[2];
sx q[2];
rz(-1.8370266) q[2];
sx q[2];
rz(0.12651786) q[2];
rz(2.8070519) q[3];
sx q[3];
rz(-2.8821475) q[3];
sx q[3];
rz(-0.87219605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8406521) q[0];
sx q[0];
rz(-2.2569188) q[0];
sx q[0];
rz(0.51710039) q[0];
rz(-3.0866947) q[1];
sx q[1];
rz(-1.6322735) q[1];
sx q[1];
rz(0.089769207) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9835998) q[0];
sx q[0];
rz(-1.9741892) q[0];
sx q[0];
rz(2.1703815) q[0];
rz(-pi) q[1];
rz(1.0099118) q[2];
sx q[2];
rz(-1.623962) q[2];
sx q[2];
rz(1.1720942) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0193737) q[1];
sx q[1];
rz(-1.8226588) q[1];
sx q[1];
rz(-1.7346965) q[1];
rz(-pi) q[2];
rz(-2.1193954) q[3];
sx q[3];
rz(-2.7226603) q[3];
sx q[3];
rz(2.0868807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6101997) q[2];
sx q[2];
rz(-2.0072082) q[2];
sx q[2];
rz(0.72719491) q[2];
rz(0.7555035) q[3];
sx q[3];
rz(-2.754039) q[3];
sx q[3];
rz(-2.9288647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1995734) q[0];
sx q[0];
rz(-1.5409536) q[0];
sx q[0];
rz(1.090747) q[0];
rz(1.749281) q[1];
sx q[1];
rz(-1.5131469) q[1];
sx q[1];
rz(-1.6319235) q[1];
rz(0.86633273) q[2];
sx q[2];
rz(-1.6523747) q[2];
sx q[2];
rz(-2.0518641) q[2];
rz(2.5531664) q[3];
sx q[3];
rz(-2.1156433) q[3];
sx q[3];
rz(2.0414203) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
