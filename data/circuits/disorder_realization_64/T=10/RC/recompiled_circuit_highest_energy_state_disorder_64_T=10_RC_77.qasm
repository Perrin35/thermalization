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
rz(-1.0567226) q[0];
sx q[0];
rz(-1.280008) q[0];
sx q[0];
rz(-0.7260538) q[0];
rz(-2.5441406) q[1];
sx q[1];
rz(-0.51550454) q[1];
sx q[1];
rz(-1.0322303) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55695483) q[0];
sx q[0];
rz(-2.224486) q[0];
sx q[0];
rz(-1.6980706) q[0];
rz(-0.039041877) q[2];
sx q[2];
rz(-1.0532951) q[2];
sx q[2];
rz(-1.3618748) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.5121442) q[1];
sx q[1];
rz(-1.8820861) q[1];
sx q[1];
rz(-1.6988847) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0809523) q[3];
sx q[3];
rz(-2.3054696) q[3];
sx q[3];
rz(1.6250798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7008179) q[2];
sx q[2];
rz(-0.60443193) q[2];
sx q[2];
rz(-3.0536998) q[2];
rz(-1.4676189) q[3];
sx q[3];
rz(-1.2940977) q[3];
sx q[3];
rz(-2.1427593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0198332) q[0];
sx q[0];
rz(-1.3061981) q[0];
sx q[0];
rz(1.4922967) q[0];
rz(2.3536033) q[1];
sx q[1];
rz(-1.183527) q[1];
sx q[1];
rz(0.96510395) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0168646) q[0];
sx q[0];
rz(-1.4788606) q[0];
sx q[0];
rz(-0.73081907) q[0];
rz(-pi) q[1];
rz(1.713322) q[2];
sx q[2];
rz(-2.3025142) q[2];
sx q[2];
rz(2.1676262) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.75747672) q[1];
sx q[1];
rz(-2.6429477) q[1];
sx q[1];
rz(1.1747141) q[1];
x q[2];
rz(-2.6718465) q[3];
sx q[3];
rz(-1.6516017) q[3];
sx q[3];
rz(1.684379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.013177055) q[2];
sx q[2];
rz(-0.33856496) q[2];
sx q[2];
rz(1.817912) q[2];
rz(2.2549818) q[3];
sx q[3];
rz(-1.9804136) q[3];
sx q[3];
rz(0.19865856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(0.24432467) q[0];
sx q[0];
rz(-0.97049814) q[0];
sx q[0];
rz(-0.6955198) q[0];
rz(2.3059402) q[1];
sx q[1];
rz(-1.4202159) q[1];
sx q[1];
rz(0.64170352) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013896522) q[0];
sx q[0];
rz(-2.3900716) q[0];
sx q[0];
rz(-1.3579382) q[0];
rz(-pi) q[1];
rz(-0.71126781) q[2];
sx q[2];
rz(-1.6984273) q[2];
sx q[2];
rz(0.39943275) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11444387) q[1];
sx q[1];
rz(-0.58747298) q[1];
sx q[1];
rz(-0.97987608) q[1];
x q[2];
rz(-0.66777457) q[3];
sx q[3];
rz(-0.80825704) q[3];
sx q[3];
rz(-0.86195743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16331638) q[2];
sx q[2];
rz(-2.0127608) q[2];
sx q[2];
rz(-1.3345435) q[2];
rz(-1.7415907) q[3];
sx q[3];
rz(-1.497523) q[3];
sx q[3];
rz(-1.1388206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9982933) q[0];
sx q[0];
rz(-2.2704953) q[0];
sx q[0];
rz(0.77348462) q[0];
rz(-0.086291226) q[1];
sx q[1];
rz(-2.6373865) q[1];
sx q[1];
rz(2.6533244) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3873169) q[0];
sx q[0];
rz(-1.5026706) q[0];
sx q[0];
rz(2.678837) q[0];
x q[1];
rz(-1.0268698) q[2];
sx q[2];
rz(-2.1928117) q[2];
sx q[2];
rz(1.1116456) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.5544071) q[1];
sx q[1];
rz(-0.76495612) q[1];
sx q[1];
rz(1.8991963) q[1];
rz(1.1639488) q[3];
sx q[3];
rz(-1.9449807) q[3];
sx q[3];
rz(-2.0896623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6733072) q[2];
sx q[2];
rz(-1.1136592) q[2];
sx q[2];
rz(0.16806531) q[2];
rz(0.42260653) q[3];
sx q[3];
rz(-0.97413078) q[3];
sx q[3];
rz(-0.15338038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76822686) q[0];
sx q[0];
rz(-0.90665561) q[0];
sx q[0];
rz(-1.9510829) q[0];
rz(-2.6137784) q[1];
sx q[1];
rz(-1.0820791) q[1];
sx q[1];
rz(-3.1324918) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21484387) q[0];
sx q[0];
rz(-2.0159449) q[0];
sx q[0];
rz(0.67589228) q[0];
x q[1];
rz(-2.3528655) q[2];
sx q[2];
rz(-0.085914748) q[2];
sx q[2];
rz(-2.9979561) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2796265) q[1];
sx q[1];
rz(-1.990296) q[1];
sx q[1];
rz(2.2395833) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9420932) q[3];
sx q[3];
rz(-1.9607984) q[3];
sx q[3];
rz(2.3671234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3147543) q[2];
sx q[2];
rz(-0.49586168) q[2];
sx q[2];
rz(-0.66519386) q[2];
rz(-2.1151755) q[3];
sx q[3];
rz(-2.0668991) q[3];
sx q[3];
rz(-0.78072602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.6219532) q[0];
sx q[0];
rz(-2.6060947) q[0];
sx q[0];
rz(-0.39805472) q[0];
rz(-1.6710501) q[1];
sx q[1];
rz(-2.2649951) q[1];
sx q[1];
rz(-2.3614531) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1127284) q[0];
sx q[0];
rz(-0.60883373) q[0];
sx q[0];
rz(-0.91797773) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8988437) q[2];
sx q[2];
rz(-1.7673337) q[2];
sx q[2];
rz(2.8145973) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8079905) q[1];
sx q[1];
rz(-1.2871075) q[1];
sx q[1];
rz(0.70244838) q[1];
rz(0.8858812) q[3];
sx q[3];
rz(-1.8198593) q[3];
sx q[3];
rz(2.0307756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2842497) q[2];
sx q[2];
rz(-2.0390859) q[2];
sx q[2];
rz(0.16172376) q[2];
rz(3.0297847) q[3];
sx q[3];
rz(-0.72607741) q[3];
sx q[3];
rz(-2.4017754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9991456) q[0];
sx q[0];
rz(-1.4302) q[0];
sx q[0];
rz(-0.69851843) q[0];
rz(2.0493719) q[1];
sx q[1];
rz(-0.75463808) q[1];
sx q[1];
rz(0.05365595) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.221401) q[0];
sx q[0];
rz(-0.86556731) q[0];
sx q[0];
rz(-0.34366519) q[0];
rz(-0.8648134) q[2];
sx q[2];
rz(-2.5891782) q[2];
sx q[2];
rz(-0.68589003) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.21660575) q[1];
sx q[1];
rz(-2.0408071) q[1];
sx q[1];
rz(-1.5108677) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7683623) q[3];
sx q[3];
rz(-0.65614349) q[3];
sx q[3];
rz(-0.11696091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1412389) q[2];
sx q[2];
rz(-0.69811368) q[2];
sx q[2];
rz(-2.9290283) q[2];
rz(0.32771787) q[3];
sx q[3];
rz(-2.1543584) q[3];
sx q[3];
rz(-1.7879965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.89606365) q[0];
sx q[0];
rz(-2.0241757) q[0];
sx q[0];
rz(-0.32972202) q[0];
rz(-0.7715191) q[1];
sx q[1];
rz(-0.78116575) q[1];
sx q[1];
rz(2.3158997) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2742608) q[0];
sx q[0];
rz(-1.4837588) q[0];
sx q[0];
rz(0.66477338) q[0];
rz(-0.97445455) q[2];
sx q[2];
rz(-1.6357517) q[2];
sx q[2];
rz(0.37298733) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4976501) q[1];
sx q[1];
rz(-1.4621328) q[1];
sx q[1];
rz(-0.76246271) q[1];
rz(1.3729893) q[3];
sx q[3];
rz(-1.8231492) q[3];
sx q[3];
rz(-2.7276602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.19994911) q[2];
sx q[2];
rz(-1.7295094) q[2];
sx q[2];
rz(1.4818209) q[2];
rz(-3.0485349) q[3];
sx q[3];
rz(-1.2991644) q[3];
sx q[3];
rz(-0.53946462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2411497) q[0];
sx q[0];
rz(-2.7261782) q[0];
sx q[0];
rz(-2.7442617) q[0];
rz(-0.98995248) q[1];
sx q[1];
rz(-1.3497458) q[1];
sx q[1];
rz(2.1053402) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66175211) q[0];
sx q[0];
rz(-1.9755198) q[0];
sx q[0];
rz(-1.5581801) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7781899) q[2];
sx q[2];
rz(-0.37523233) q[2];
sx q[2];
rz(2.4979532) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.81779382) q[1];
sx q[1];
rz(-1.7897871) q[1];
sx q[1];
rz(1.4794255) q[1];
rz(2.1734222) q[3];
sx q[3];
rz(-2.535248) q[3];
sx q[3];
rz(-1.2392288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.022126023) q[2];
sx q[2];
rz(-0.72212044) q[2];
sx q[2];
rz(1.3163346) q[2];
rz(0.25257603) q[3];
sx q[3];
rz(-1.794869) q[3];
sx q[3];
rz(-0.36309567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0040358) q[0];
sx q[0];
rz(-0.79020774) q[0];
sx q[0];
rz(-2.9055415) q[0];
rz(3.0421742) q[1];
sx q[1];
rz(-0.75629083) q[1];
sx q[1];
rz(-1.258446) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6981068) q[0];
sx q[0];
rz(-0.56523609) q[0];
sx q[0];
rz(-0.14552571) q[0];
rz(-0.66345556) q[2];
sx q[2];
rz(-1.1642312) q[2];
sx q[2];
rz(0.50944177) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.26791456) q[1];
sx q[1];
rz(-2.2482052) q[1];
sx q[1];
rz(1.4235086) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80983617) q[3];
sx q[3];
rz(-2.5812529) q[3];
sx q[3];
rz(2.5446825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1666169) q[2];
sx q[2];
rz(-0.48305837) q[2];
sx q[2];
rz(2.100259) q[2];
rz(-2.5552022) q[3];
sx q[3];
rz(-0.47724637) q[3];
sx q[3];
rz(-2.3311116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36622421) q[0];
sx q[0];
rz(-2.0851705) q[0];
sx q[0];
rz(2.0023517) q[0];
rz(2.1495023) q[1];
sx q[1];
rz(-1.5403668) q[1];
sx q[1];
rz(-1.5425727) q[1];
rz(-0.31789657) q[2];
sx q[2];
rz(-1.7902725) q[2];
sx q[2];
rz(-0.27223311) q[2];
rz(0.6057779) q[3];
sx q[3];
rz(-0.68690261) q[3];
sx q[3];
rz(1.5920832) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
