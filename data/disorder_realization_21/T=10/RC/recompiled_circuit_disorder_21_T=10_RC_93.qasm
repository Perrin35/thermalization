OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.055846) q[0];
sx q[0];
rz(6.364967) q[0];
sx q[0];
rz(9.9262417) q[0];
rz(-1.6429098) q[1];
sx q[1];
rz(-0.39615762) q[1];
sx q[1];
rz(0.3224386) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6591588) q[0];
sx q[0];
rz(-1.4855488) q[0];
sx q[0];
rz(-1.2486588) q[0];
rz(0.88959496) q[2];
sx q[2];
rz(-2.0686364) q[2];
sx q[2];
rz(-1.7878469) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8391708) q[1];
sx q[1];
rz(-1.2245373) q[1];
sx q[1];
rz(-2.195921) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93482165) q[3];
sx q[3];
rz(-1.9536195) q[3];
sx q[3];
rz(-2.8179907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.50513187) q[2];
sx q[2];
rz(-0.59288609) q[2];
sx q[2];
rz(2.5855529) q[2];
rz(0.83267823) q[3];
sx q[3];
rz(-1.6502389) q[3];
sx q[3];
rz(-2.1957943) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44822025) q[0];
sx q[0];
rz(-1.4602666) q[0];
sx q[0];
rz(-2.9843176) q[0];
rz(-0.26113025) q[1];
sx q[1];
rz(-1.3477247) q[1];
sx q[1];
rz(-0.10903407) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029023829) q[0];
sx q[0];
rz(-1.3431664) q[0];
sx q[0];
rz(-2.9665222) q[0];
rz(-0.26308665) q[2];
sx q[2];
rz(-1.7619942) q[2];
sx q[2];
rz(2.4441602) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.82169689) q[1];
sx q[1];
rz(-0.79274717) q[1];
sx q[1];
rz(2.4577623) q[1];
rz(-pi) q[2];
rz(-0.29626366) q[3];
sx q[3];
rz(-0.54013541) q[3];
sx q[3];
rz(1.4512645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9033501) q[2];
sx q[2];
rz(-1.976333) q[2];
sx q[2];
rz(-1.8781352) q[2];
rz(2.8144828) q[3];
sx q[3];
rz(-1.5644904) q[3];
sx q[3];
rz(1.2143149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-3.0771714) q[0];
sx q[0];
rz(-0.049296878) q[0];
sx q[0];
rz(1.7984614) q[0];
rz(2.893977) q[1];
sx q[1];
rz(-2.394948) q[1];
sx q[1];
rz(-0.48167357) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7307229) q[0];
sx q[0];
rz(-2.0559089) q[0];
sx q[0];
rz(-2.7184125) q[0];
rz(-2.8086497) q[2];
sx q[2];
rz(-0.98368401) q[2];
sx q[2];
rz(-1.743403) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.86320089) q[1];
sx q[1];
rz(-2.26483) q[1];
sx q[1];
rz(0.73514003) q[1];
rz(-0.89110156) q[3];
sx q[3];
rz(-2.292233) q[3];
sx q[3];
rz(-0.47618714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8032916) q[2];
sx q[2];
rz(-0.81739601) q[2];
sx q[2];
rz(0.49989191) q[2];
rz(2.5806184) q[3];
sx q[3];
rz(-1.8818972) q[3];
sx q[3];
rz(-1.6803754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19701476) q[0];
sx q[0];
rz(-2.975583) q[0];
sx q[0];
rz(-0.552185) q[0];
rz(-1.588297) q[1];
sx q[1];
rz(-2.242656) q[1];
sx q[1];
rz(1.8968556) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23918505) q[0];
sx q[0];
rz(-2.8021078) q[0];
sx q[0];
rz(3.1367338) q[0];
rz(3.0110637) q[2];
sx q[2];
rz(-1.7703238) q[2];
sx q[2];
rz(2.5752657) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1340449) q[1];
sx q[1];
rz(-1.3589077) q[1];
sx q[1];
rz(-2.120554) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57085412) q[3];
sx q[3];
rz(-1.7613162) q[3];
sx q[3];
rz(-1.1754456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2923979) q[2];
sx q[2];
rz(-1.8820102) q[2];
sx q[2];
rz(-1.9909031) q[2];
rz(-1.6644647) q[3];
sx q[3];
rz(-1.632558) q[3];
sx q[3];
rz(-0.48294827) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0634336) q[0];
sx q[0];
rz(-0.76197356) q[0];
sx q[0];
rz(-0.081469014) q[0];
rz(3.07913) q[1];
sx q[1];
rz(-1.1413347) q[1];
sx q[1];
rz(1.5030456) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7446639) q[0];
sx q[0];
rz(-1.6984807) q[0];
sx q[0];
rz(0.98608195) q[0];
rz(-pi) q[1];
rz(1.6790752) q[2];
sx q[2];
rz(-1.6533972) q[2];
sx q[2];
rz(0.29512197) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.45080966) q[1];
sx q[1];
rz(-2.8125617) q[1];
sx q[1];
rz(-1.8915218) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53346975) q[3];
sx q[3];
rz(-0.57538486) q[3];
sx q[3];
rz(-1.7076147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6333255) q[2];
sx q[2];
rz(-2.1990364) q[2];
sx q[2];
rz(1.903669) q[2];
rz(2.0189019) q[3];
sx q[3];
rz(-2.4653547) q[3];
sx q[3];
rz(0.52156633) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6102819) q[0];
sx q[0];
rz(-2.1370482) q[0];
sx q[0];
rz(2.8748728) q[0];
rz(-0.56089127) q[1];
sx q[1];
rz(-1.8436878) q[1];
sx q[1];
rz(2.3430603) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7245367) q[0];
sx q[0];
rz(-2.8310611) q[0];
sx q[0];
rz(1.0945555) q[0];
x q[1];
rz(2.550225) q[2];
sx q[2];
rz(-2.5632576) q[2];
sx q[2];
rz(0.27507281) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.861607) q[1];
sx q[1];
rz(-1.0000739) q[1];
sx q[1];
rz(0.062203783) q[1];
rz(0.60243209) q[3];
sx q[3];
rz(-0.89655399) q[3];
sx q[3];
rz(2.7979421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.16053998) q[2];
sx q[2];
rz(-1.2489698) q[2];
sx q[2];
rz(2.7748761) q[2];
rz(1.8803053) q[3];
sx q[3];
rz(-1.6882608) q[3];
sx q[3];
rz(3.0453851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7325608) q[0];
sx q[0];
rz(-0.92027396) q[0];
sx q[0];
rz(-0.60638705) q[0];
rz(-0.19730332) q[1];
sx q[1];
rz(-1.1261255) q[1];
sx q[1];
rz(0.46404776) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2460829) q[0];
sx q[0];
rz(-1.1469139) q[0];
sx q[0];
rz(-2.258582) q[0];
rz(-pi) q[1];
rz(0.020178528) q[2];
sx q[2];
rz(-1.8115461) q[2];
sx q[2];
rz(-1.9764331) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6021767) q[1];
sx q[1];
rz(-0.5836986) q[1];
sx q[1];
rz(-2.1392439) q[1];
rz(0.32902284) q[3];
sx q[3];
rz(-0.24917069) q[3];
sx q[3];
rz(-0.84013018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.93418926) q[2];
sx q[2];
rz(-2.1384017) q[2];
sx q[2];
rz(0.25804538) q[2];
rz(1.1856273) q[3];
sx q[3];
rz(-1.6262755) q[3];
sx q[3];
rz(0.0035088249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.946452) q[0];
sx q[0];
rz(-1.2807245) q[0];
sx q[0];
rz(-0.38129693) q[0];
rz(-0.095245846) q[1];
sx q[1];
rz(-2.1691599) q[1];
sx q[1];
rz(1.4415178) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9593175) q[0];
sx q[0];
rz(-1.5858985) q[0];
sx q[0];
rz(-3.0781151) q[0];
rz(2.5289815) q[2];
sx q[2];
rz(-1.6022041) q[2];
sx q[2];
rz(-2.999246) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0563593) q[1];
sx q[1];
rz(-1.4409522) q[1];
sx q[1];
rz(-1.4443881) q[1];
rz(-1.6137245) q[3];
sx q[3];
rz(-1.9044442) q[3];
sx q[3];
rz(-1.4294525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9986481) q[2];
sx q[2];
rz(-0.412985) q[2];
sx q[2];
rz(-2.9150035) q[2];
rz(-0.44858027) q[3];
sx q[3];
rz(-1.5357163) q[3];
sx q[3];
rz(-2.3118238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7609693) q[0];
sx q[0];
rz(-0.7614823) q[0];
sx q[0];
rz(-1.7425849) q[0];
rz(0.31708583) q[1];
sx q[1];
rz(-1.6665019) q[1];
sx q[1];
rz(0.98659602) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50842677) q[0];
sx q[0];
rz(-1.0796483) q[0];
sx q[0];
rz(1.7500061) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87395845) q[2];
sx q[2];
rz(-1.9947589) q[2];
sx q[2];
rz(1.6800113) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7762737) q[1];
sx q[1];
rz(-2.1734108) q[1];
sx q[1];
rz(-1.7425294) q[1];
rz(-pi) q[2];
rz(1.1986198) q[3];
sx q[3];
rz(-0.44789207) q[3];
sx q[3];
rz(-0.4263634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2150779) q[2];
sx q[2];
rz(-2.4145917) q[2];
sx q[2];
rz(-2.731936) q[2];
rz(0.26327291) q[3];
sx q[3];
rz(-1.3116838) q[3];
sx q[3];
rz(-0.51945654) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7423994) q[0];
sx q[0];
rz(-0.078646794) q[0];
sx q[0];
rz(1.7364527) q[0];
rz(-2.3204904) q[1];
sx q[1];
rz(-2.2228873) q[1];
sx q[1];
rz(-1.7260889) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4840568) q[0];
sx q[0];
rz(-2.7298379) q[0];
sx q[0];
rz(0.62396892) q[0];
x q[1];
rz(1.2070933) q[2];
sx q[2];
rz(-1.8038245) q[2];
sx q[2];
rz(-1.2365014) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0714598) q[1];
sx q[1];
rz(-1.7839285) q[1];
sx q[1];
rz(0.15503426) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0584162) q[3];
sx q[3];
rz(-1.7184162) q[3];
sx q[3];
rz(0.83422134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4225509) q[2];
sx q[2];
rz(-2.8406403) q[2];
sx q[2];
rz(0.12410513) q[2];
rz(0.96578807) q[3];
sx q[3];
rz(-1.4744759) q[3];
sx q[3];
rz(2.4805099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60349764) q[0];
sx q[0];
rz(-0.24833831) q[0];
sx q[0];
rz(-0.86059358) q[0];
rz(-2.8339236) q[1];
sx q[1];
rz(-1.888702) q[1];
sx q[1];
rz(-1.9370334) q[1];
rz(2.5406191) q[2];
sx q[2];
rz(-2.8291611) q[2];
sx q[2];
rz(-0.57401382) q[2];
rz(-1.4205167) q[3];
sx q[3];
rz(-1.0192623) q[3];
sx q[3];
rz(-2.7912959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];