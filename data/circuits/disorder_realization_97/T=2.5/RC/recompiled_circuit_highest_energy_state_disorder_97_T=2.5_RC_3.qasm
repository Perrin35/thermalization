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
rz(-0.66145575) q[0];
sx q[0];
rz(-0.28132004) q[0];
sx q[0];
rz(1.1433648) q[0];
rz(0.65302628) q[1];
sx q[1];
rz(4.9622494) q[1];
sx q[1];
rz(9.4402037) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.832725) q[0];
sx q[0];
rz(-1.0842396) q[0];
sx q[0];
rz(2.1910446) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77818971) q[2];
sx q[2];
rz(-2.3986849) q[2];
sx q[2];
rz(-1.6011432) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.28492883) q[1];
sx q[1];
rz(-0.94051111) q[1];
sx q[1];
rz(0.50215118) q[1];
x q[2];
rz(-1.3728834) q[3];
sx q[3];
rz(-2.1002401) q[3];
sx q[3];
rz(2.6967665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.780484) q[2];
sx q[2];
rz(-0.0094272308) q[2];
sx q[2];
rz(2.0375605) q[2];
rz(2.5822254) q[3];
sx q[3];
rz(-2.1909824) q[3];
sx q[3];
rz(-0.88421384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1169432) q[0];
sx q[0];
rz(-0.58498061) q[0];
sx q[0];
rz(0.90644932) q[0];
rz(-2.8001884) q[1];
sx q[1];
rz(-2.2661792) q[1];
sx q[1];
rz(2.7934449) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6432836) q[0];
sx q[0];
rz(-1.8178682) q[0];
sx q[0];
rz(-0.31881316) q[0];
rz(-pi) q[1];
rz(-2.1402713) q[2];
sx q[2];
rz(-3.0932326) q[2];
sx q[2];
rz(2.0404301) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9056919) q[1];
sx q[1];
rz(-2.5254619) q[1];
sx q[1];
rz(1.4409164) q[1];
rz(0.23032974) q[3];
sx q[3];
rz(-0.84127142) q[3];
sx q[3];
rz(-1.426633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.97588313) q[2];
sx q[2];
rz(-0.3419064) q[2];
sx q[2];
rz(-3.0561786) q[2];
rz(0.092770569) q[3];
sx q[3];
rz(-2.2723891) q[3];
sx q[3];
rz(2.2658277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77466011) q[0];
sx q[0];
rz(-2.7800738) q[0];
sx q[0];
rz(2.7969978) q[0];
rz(1.5917646) q[1];
sx q[1];
rz(-1.7235618) q[1];
sx q[1];
rz(1.4580911) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34209308) q[0];
sx q[0];
rz(-1.5342426) q[0];
sx q[0];
rz(0.7176368) q[0];
rz(2.0877203) q[2];
sx q[2];
rz(-1.6224684) q[2];
sx q[2];
rz(-1.8926562) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5418936) q[1];
sx q[1];
rz(-2.8435282) q[1];
sx q[1];
rz(0.96082176) q[1];
rz(-pi) q[2];
rz(-0.97403861) q[3];
sx q[3];
rz(-2.0152053) q[3];
sx q[3];
rz(-0.35154464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6452667) q[2];
sx q[2];
rz(-0.93924773) q[2];
sx q[2];
rz(-1.1305031) q[2];
rz(0.94275236) q[3];
sx q[3];
rz(-1.0470942) q[3];
sx q[3];
rz(-1.2070791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4261674) q[0];
sx q[0];
rz(-1.796145) q[0];
sx q[0];
rz(-1.6834393) q[0];
rz(-1.0484877) q[1];
sx q[1];
rz(-1.6494992) q[1];
sx q[1];
rz(2.6715211) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9662358) q[0];
sx q[0];
rz(-0.58307532) q[0];
sx q[0];
rz(-2.6158041) q[0];
rz(-pi) q[1];
rz(2.0087035) q[2];
sx q[2];
rz(-0.57269579) q[2];
sx q[2];
rz(-1.1066135) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.64420358) q[1];
sx q[1];
rz(-2.4371689) q[1];
sx q[1];
rz(0.299244) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80959971) q[3];
sx q[3];
rz(-2.2815846) q[3];
sx q[3];
rz(-0.87464911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5423535) q[2];
sx q[2];
rz(-0.75672823) q[2];
sx q[2];
rz(-2.2139464) q[2];
rz(1.2841691) q[3];
sx q[3];
rz(-1.410306) q[3];
sx q[3];
rz(-0.94902432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1272142) q[0];
sx q[0];
rz(-0.40507409) q[0];
sx q[0];
rz(0.48144427) q[0];
rz(-0.97775793) q[1];
sx q[1];
rz(-0.6441741) q[1];
sx q[1];
rz(1.0747304) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6310997) q[0];
sx q[0];
rz(-1.8429776) q[0];
sx q[0];
rz(-2.0049606) q[0];
rz(-1.85249) q[2];
sx q[2];
rz(-1.2810724) q[2];
sx q[2];
rz(0.85838028) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2664484) q[1];
sx q[1];
rz(-1.9423264) q[1];
sx q[1];
rz(0.2569916) q[1];
rz(-pi) q[2];
rz(-0.67258622) q[3];
sx q[3];
rz(-2.6574316) q[3];
sx q[3];
rz(-1.3369833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1017199) q[2];
sx q[2];
rz(-1.1168672) q[2];
sx q[2];
rz(2.5679585) q[2];
rz(-1.4839858) q[3];
sx q[3];
rz(-2.0935757) q[3];
sx q[3];
rz(-0.60605961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1178591) q[0];
sx q[0];
rz(-0.25256279) q[0];
sx q[0];
rz(0.31164393) q[0];
rz(0.32870865) q[1];
sx q[1];
rz(-1.8054104) q[1];
sx q[1];
rz(0.74105826) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5609076) q[0];
sx q[0];
rz(-1.1162762) q[0];
sx q[0];
rz(0.88991388) q[0];
rz(-pi) q[1];
rz(-1.7996721) q[2];
sx q[2];
rz(-1.5107578) q[2];
sx q[2];
rz(-2.6031074) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0892798) q[1];
sx q[1];
rz(-1.5075753) q[1];
sx q[1];
rz(0.88892816) q[1];
rz(-0.93993069) q[3];
sx q[3];
rz(-1.5700911) q[3];
sx q[3];
rz(3.0169064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7753503) q[2];
sx q[2];
rz(-0.58997184) q[2];
sx q[2];
rz(-0.74014202) q[2];
rz(0.38672334) q[3];
sx q[3];
rz(-0.72606483) q[3];
sx q[3];
rz(2.6630785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9750403) q[0];
sx q[0];
rz(-0.98169011) q[0];
sx q[0];
rz(2.8983086) q[0];
rz(2.2443306) q[1];
sx q[1];
rz(-1.5586531) q[1];
sx q[1];
rz(2.3975587) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76975856) q[0];
sx q[0];
rz(-2.7427865) q[0];
sx q[0];
rz(-0.21488667) q[0];
rz(-pi) q[1];
rz(0.93397385) q[2];
sx q[2];
rz(-1.8806071) q[2];
sx q[2];
rz(-1.241809) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1330382) q[1];
sx q[1];
rz(-2.4835476) q[1];
sx q[1];
rz(-2.1997019) q[1];
rz(-pi) q[2];
rz(-0.18455581) q[3];
sx q[3];
rz(-2.089644) q[3];
sx q[3];
rz(2.542582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9545142) q[2];
sx q[2];
rz(-2.3176471) q[2];
sx q[2];
rz(0.0079060923) q[2];
rz(-1.4963957) q[3];
sx q[3];
rz(-0.073286101) q[3];
sx q[3];
rz(2.6325398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.024648333) q[0];
sx q[0];
rz(-0.032289676) q[0];
sx q[0];
rz(0.13667983) q[0];
rz(1.3487123) q[1];
sx q[1];
rz(-1.8486479) q[1];
sx q[1];
rz(-0.50833702) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0011372066) q[0];
sx q[0];
rz(-1.1759725) q[0];
sx q[0];
rz(3.0009901) q[0];
rz(-pi) q[1];
rz(1.2042785) q[2];
sx q[2];
rz(-2.4319639) q[2];
sx q[2];
rz(0.21096551) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6189948) q[1];
sx q[1];
rz(-0.77678372) q[1];
sx q[1];
rz(-2.0911699) q[1];
x q[2];
rz(1.6898722) q[3];
sx q[3];
rz(-1.649646) q[3];
sx q[3];
rz(-0.54311968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6792182) q[2];
sx q[2];
rz(-2.4269673) q[2];
sx q[2];
rz(-0.057961658) q[2];
rz(-2.1884749) q[3];
sx q[3];
rz(-1.0473731) q[3];
sx q[3];
rz(-2.5061149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6111095) q[0];
sx q[0];
rz(-1.4169175) q[0];
sx q[0];
rz(0.86607754) q[0];
rz(-1.3653612) q[1];
sx q[1];
rz(-1.5581286) q[1];
sx q[1];
rz(2.3376215) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83528501) q[0];
sx q[0];
rz(-0.87125766) q[0];
sx q[0];
rz(0.37658306) q[0];
rz(1.1273659) q[2];
sx q[2];
rz(-0.83258234) q[2];
sx q[2];
rz(-1.2633367) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.44195804) q[1];
sx q[1];
rz(-0.82773877) q[1];
sx q[1];
rz(2.1585967) q[1];
rz(1.9744121) q[3];
sx q[3];
rz(-1.5506553) q[3];
sx q[3];
rz(-0.27475629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0552401) q[2];
sx q[2];
rz(-0.1487727) q[2];
sx q[2];
rz(-2.0894334) q[2];
rz(2.2822028) q[3];
sx q[3];
rz(-0.79901564) q[3];
sx q[3];
rz(0.94256443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6732442) q[0];
sx q[0];
rz(-0.66272658) q[0];
sx q[0];
rz(-2.7224139) q[0];
rz(-3.1255417) q[1];
sx q[1];
rz(-1.586986) q[1];
sx q[1];
rz(3.0174461) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9036839) q[0];
sx q[0];
rz(-2.5669837) q[0];
sx q[0];
rz(1.8540856) q[0];
rz(-pi) q[1];
rz(2.2474278) q[2];
sx q[2];
rz(-2.8196206) q[2];
sx q[2];
rz(1.8338667) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.649227) q[1];
sx q[1];
rz(-1.320663) q[1];
sx q[1];
rz(1.753128) q[1];
rz(-0.18245635) q[3];
sx q[3];
rz(-2.2130551) q[3];
sx q[3];
rz(0.69891847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.19196709) q[2];
sx q[2];
rz(-0.73835915) q[2];
sx q[2];
rz(0.16853608) q[2];
rz(-1.4847697) q[3];
sx q[3];
rz(-0.46983457) q[3];
sx q[3];
rz(2.7114939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94854245) q[0];
sx q[0];
rz(-0.47946231) q[0];
sx q[0];
rz(-0.23647501) q[0];
rz(-2.3169658) q[1];
sx q[1];
rz(-1.6025447) q[1];
sx q[1];
rz(-1.2784169) q[1];
rz(2.3659351) q[2];
sx q[2];
rz(-2.0022598) q[2];
sx q[2];
rz(-0.74874292) q[2];
rz(-1.6705728) q[3];
sx q[3];
rz(-2.592516) q[3];
sx q[3];
rz(1.9883131) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
