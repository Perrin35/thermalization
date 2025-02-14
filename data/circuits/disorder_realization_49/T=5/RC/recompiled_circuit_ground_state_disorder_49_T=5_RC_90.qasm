OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9743118) q[0];
sx q[0];
rz(-0.33742961) q[0];
sx q[0];
rz(0.20198241) q[0];
rz(2.1972411) q[1];
sx q[1];
rz(-0.75466067) q[1];
sx q[1];
rz(1.4442297) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9044735) q[0];
sx q[0];
rz(-0.54330641) q[0];
sx q[0];
rz(-2.7271557) q[0];
x q[1];
rz(-0.10798048) q[2];
sx q[2];
rz(-1.3657346) q[2];
sx q[2];
rz(2.972796) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3814427) q[1];
sx q[1];
rz(-2.0408148) q[1];
sx q[1];
rz(-0.85531855) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.045363868) q[3];
sx q[3];
rz(-1.9352018) q[3];
sx q[3];
rz(1.7045329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.28295383) q[2];
sx q[2];
rz(-0.084246548) q[2];
sx q[2];
rz(1.5343182) q[2];
rz(-2.9547847) q[3];
sx q[3];
rz(-0.45972937) q[3];
sx q[3];
rz(1.648858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9635791) q[0];
sx q[0];
rz(-2.3602965) q[0];
sx q[0];
rz(2.436893) q[0];
rz(1.0164545) q[1];
sx q[1];
rz(-1.9489138) q[1];
sx q[1];
rz(1.1526398) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64166028) q[0];
sx q[0];
rz(-1.0295273) q[0];
sx q[0];
rz(2.9425147) q[0];
rz(-2.1179885) q[2];
sx q[2];
rz(-2.0788801) q[2];
sx q[2];
rz(-1.2244727) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1981467) q[1];
sx q[1];
rz(-2.0753161) q[1];
sx q[1];
rz(-2.2322502) q[1];
rz(-0.74781466) q[3];
sx q[3];
rz(-1.4749817) q[3];
sx q[3];
rz(0.15296061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.54027259) q[2];
sx q[2];
rz(-2.2728964) q[2];
sx q[2];
rz(1.3199838) q[2];
rz(-0.071062239) q[3];
sx q[3];
rz(-0.080852121) q[3];
sx q[3];
rz(-2.0474056) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0367301) q[0];
sx q[0];
rz(-2.7920089) q[0];
sx q[0];
rz(-2.2678251) q[0];
rz(1.8050487) q[1];
sx q[1];
rz(-2.4599894) q[1];
sx q[1];
rz(-2.2866586) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0655496) q[0];
sx q[0];
rz(-2.4618755) q[0];
sx q[0];
rz(1.9125841) q[0];
x q[1];
rz(-2.1941691) q[2];
sx q[2];
rz(-0.50259841) q[2];
sx q[2];
rz(-1.1029152) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1347772) q[1];
sx q[1];
rz(-1.8726649) q[1];
sx q[1];
rz(1.6848438) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7822262) q[3];
sx q[3];
rz(-0.99726652) q[3];
sx q[3];
rz(1.9952967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9285589) q[2];
sx q[2];
rz(-1.1308257) q[2];
sx q[2];
rz(1.0553168) q[2];
rz(2.1999551) q[3];
sx q[3];
rz(-2.0703273) q[3];
sx q[3];
rz(-2.1987703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68614352) q[0];
sx q[0];
rz(-2.2930155) q[0];
sx q[0];
rz(-0.69547478) q[0];
rz(1.4675568) q[1];
sx q[1];
rz(-1.2662042) q[1];
sx q[1];
rz(-3.0470336) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5887711) q[0];
sx q[0];
rz(-1.2929734) q[0];
sx q[0];
rz(0.050006067) q[0];
x q[1];
rz(-0.84708656) q[2];
sx q[2];
rz(-2.0485592) q[2];
sx q[2];
rz(1.0003045) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.70909158) q[1];
sx q[1];
rz(-1.7686262) q[1];
sx q[1];
rz(1.2835931) q[1];
x q[2];
rz(-1.2572722) q[3];
sx q[3];
rz(-0.71259585) q[3];
sx q[3];
rz(-0.37171504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.62119421) q[2];
sx q[2];
rz(-0.79402557) q[2];
sx q[2];
rz(0.76148477) q[2];
rz(-0.17101184) q[3];
sx q[3];
rz(-1.2820425) q[3];
sx q[3];
rz(3.0448992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3068202) q[0];
sx q[0];
rz(-0.67245317) q[0];
sx q[0];
rz(-2.6493454) q[0];
rz(-1.8222088) q[1];
sx q[1];
rz(-1.4232114) q[1];
sx q[1];
rz(-0.48670235) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20535417) q[0];
sx q[0];
rz(-1.1726609) q[0];
sx q[0];
rz(-1.7004844) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7560744) q[2];
sx q[2];
rz(-0.80759128) q[2];
sx q[2];
rz(2.453192) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4449205) q[1];
sx q[1];
rz(-1.2081523) q[1];
sx q[1];
rz(1.1128694) q[1];
x q[2];
rz(-1.9697804) q[3];
sx q[3];
rz(-2.5907142) q[3];
sx q[3];
rz(-0.18607947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1897366) q[2];
sx q[2];
rz(-0.8437914) q[2];
sx q[2];
rz(-1.0978511) q[2];
rz(2.3073933) q[3];
sx q[3];
rz(-0.67535496) q[3];
sx q[3];
rz(-1.7043017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6247691) q[0];
sx q[0];
rz(-0.56684816) q[0];
sx q[0];
rz(-3.095742) q[0];
rz(-0.74371964) q[1];
sx q[1];
rz(-2.7727978) q[1];
sx q[1];
rz(-3.0812982) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3876095) q[0];
sx q[0];
rz(-1.930325) q[0];
sx q[0];
rz(0.44937581) q[0];
rz(-pi) q[1];
rz(-0.8745114) q[2];
sx q[2];
rz(-1.2605476) q[2];
sx q[2];
rz(-0.53521148) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.454596) q[1];
sx q[1];
rz(-2.7424062) q[1];
sx q[1];
rz(2.1737425) q[1];
rz(-pi) q[2];
rz(1.4521818) q[3];
sx q[3];
rz(-2.4595692) q[3];
sx q[3];
rz(-0.31957808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.22659773) q[2];
sx q[2];
rz(-1.3264341) q[2];
sx q[2];
rz(2.162852) q[2];
rz(-1.2616875) q[3];
sx q[3];
rz(-1.0190957) q[3];
sx q[3];
rz(-2.252388) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17734811) q[0];
sx q[0];
rz(-1.526399) q[0];
sx q[0];
rz(-0.88441315) q[0];
rz(2.2824967) q[1];
sx q[1];
rz(-1.1235378) q[1];
sx q[1];
rz(-0.20800796) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1387062) q[0];
sx q[0];
rz(-2.0660668) q[0];
sx q[0];
rz(1.9751092) q[0];
x q[1];
rz(-2.538717) q[2];
sx q[2];
rz(-1.2211717) q[2];
sx q[2];
rz(2.5536551) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2629548) q[1];
sx q[1];
rz(-2.4391101) q[1];
sx q[1];
rz(1.9874057) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0741803) q[3];
sx q[3];
rz(-1.5690104) q[3];
sx q[3];
rz(-1.3942277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.38006833) q[2];
sx q[2];
rz(-0.91975776) q[2];
sx q[2];
rz(-2.5999787) q[2];
rz(1.2540865) q[3];
sx q[3];
rz(-1.5518291) q[3];
sx q[3];
rz(-0.90887535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33614531) q[0];
sx q[0];
rz(-2.8271524) q[0];
sx q[0];
rz(-2.050515) q[0];
rz(-0.82487851) q[1];
sx q[1];
rz(-2.7179317) q[1];
sx q[1];
rz(2.0687912) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7110853) q[0];
sx q[0];
rz(-0.46784624) q[0];
sx q[0];
rz(-0.43456315) q[0];
x q[1];
rz(-0.013117803) q[2];
sx q[2];
rz(-1.5481536) q[2];
sx q[2];
rz(-2.4215339) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.78847105) q[1];
sx q[1];
rz(-1.7040729) q[1];
sx q[1];
rz(-0.2443831) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0951525) q[3];
sx q[3];
rz(-0.45733157) q[3];
sx q[3];
rz(-0.90378896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.23952809) q[2];
sx q[2];
rz(-0.54486474) q[2];
sx q[2];
rz(-2.876335) q[2];
rz(-0.15051633) q[3];
sx q[3];
rz(-2.033332) q[3];
sx q[3];
rz(0.36293852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4758258) q[0];
sx q[0];
rz(-0.095223991) q[0];
sx q[0];
rz(-1.6640523) q[0];
rz(2.3382969) q[1];
sx q[1];
rz(-1.7684312) q[1];
sx q[1];
rz(1.0172179) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2870713) q[0];
sx q[0];
rz(-1.4421363) q[0];
sx q[0];
rz(1.3657938) q[0];
rz(-pi) q[1];
rz(1.3147589) q[2];
sx q[2];
rz(-2.0190329) q[2];
sx q[2];
rz(-2.2729276) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.13297824) q[1];
sx q[1];
rz(-2.6854351) q[1];
sx q[1];
rz(0.71059473) q[1];
x q[2];
rz(-0.32633304) q[3];
sx q[3];
rz(-1.2830955) q[3];
sx q[3];
rz(-0.74244754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9251755) q[2];
sx q[2];
rz(-1.5745796) q[2];
sx q[2];
rz(0.65640059) q[2];
rz(-2.9260855) q[3];
sx q[3];
rz(-1.7301205) q[3];
sx q[3];
rz(-1.4490674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086275252) q[0];
sx q[0];
rz(-1.7973987) q[0];
sx q[0];
rz(2.2917746) q[0];
rz(0.96018106) q[1];
sx q[1];
rz(-0.4147059) q[1];
sx q[1];
rz(2.2023831) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4681743) q[0];
sx q[0];
rz(-2.3395754) q[0];
sx q[0];
rz(-1.2002263) q[0];
rz(-pi) q[1];
rz(1.1562111) q[2];
sx q[2];
rz(-1.5778981) q[2];
sx q[2];
rz(-2.9455001) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.43329217) q[1];
sx q[1];
rz(-2.1850249) q[1];
sx q[1];
rz(2.0310165) q[1];
rz(-0.47795145) q[3];
sx q[3];
rz(-2.047973) q[3];
sx q[3];
rz(1.3665703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5822997) q[2];
sx q[2];
rz(-1.7377995) q[2];
sx q[2];
rz(0.55029184) q[2];
rz(-0.127921) q[3];
sx q[3];
rz(-3.0030799) q[3];
sx q[3];
rz(0.96854717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6561103) q[0];
sx q[0];
rz(-2.4714097) q[0];
sx q[0];
rz(2.4708268) q[0];
rz(2.7863964) q[1];
sx q[1];
rz(-1.4757481) q[1];
sx q[1];
rz(-0.3955985) q[1];
rz(1.7110466) q[2];
sx q[2];
rz(-1.4288908) q[2];
sx q[2];
rz(-2.6960001) q[2];
rz(1.1200503) q[3];
sx q[3];
rz(-1.7786296) q[3];
sx q[3];
rz(2.8854388) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
