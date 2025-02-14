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
rz(-0.94435159) q[1];
sx q[1];
rz(-2.386932) q[1];
sx q[1];
rz(1.6973629) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2371191) q[0];
sx q[0];
rz(-0.54330641) q[0];
sx q[0];
rz(-0.41443698) q[0];
rz(-pi) q[1];
rz(-2.0488705) q[2];
sx q[2];
rz(-0.23140027) q[2];
sx q[2];
rz(-0.3203985) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4718402) q[1];
sx q[1];
rz(-2.3089611) q[1];
sx q[1];
rz(0.91186146) q[1];
rz(-pi) q[2];
rz(-1.6891278) q[3];
sx q[3];
rz(-2.7744996) q[3];
sx q[3];
rz(1.5637507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8586388) q[2];
sx q[2];
rz(-3.0573461) q[2];
sx q[2];
rz(1.5343182) q[2];
rz(2.9547847) q[3];
sx q[3];
rz(-2.6818633) q[3];
sx q[3];
rz(1.648858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.1780136) q[0];
sx q[0];
rz(-0.78129617) q[0];
sx q[0];
rz(-2.436893) q[0];
rz(-1.0164545) q[1];
sx q[1];
rz(-1.9489138) q[1];
sx q[1];
rz(-1.1526398) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64166028) q[0];
sx q[0];
rz(-1.0295273) q[0];
sx q[0];
rz(-0.199078) q[0];
rz(-pi) q[1];
rz(0.75147475) q[2];
sx q[2];
rz(-0.72859366) q[2];
sx q[2];
rz(2.1214243) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.94344596) q[1];
sx q[1];
rz(-1.0662765) q[1];
sx q[1];
rz(2.2322502) q[1];
rz(-pi) q[2];
rz(3.0011938) q[3];
sx q[3];
rz(-2.388846) q[3];
sx q[3];
rz(-1.6209768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.54027259) q[2];
sx q[2];
rz(-2.2728964) q[2];
sx q[2];
rz(1.3199838) q[2];
rz(-3.0705304) q[3];
sx q[3];
rz(-3.0607405) q[3];
sx q[3];
rz(1.094187) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0367301) q[0];
sx q[0];
rz(-2.7920089) q[0];
sx q[0];
rz(2.2678251) q[0];
rz(1.3365439) q[1];
sx q[1];
rz(-0.68160325) q[1];
sx q[1];
rz(0.85493404) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6470454) q[0];
sx q[0];
rz(-2.2046007) q[0];
sx q[0];
rz(0.26453544) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31051113) q[2];
sx q[2];
rz(-1.9726255) q[2];
sx q[2];
rz(-2.7257811) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3745649) q[1];
sx q[1];
rz(-0.32207707) q[1];
sx q[1];
rz(2.7912071) q[1];
rz(-pi) q[2];
rz(-1.072238) q[3];
sx q[3];
rz(-0.66594687) q[3];
sx q[3];
rz(1.7518821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9285589) q[2];
sx q[2];
rz(-2.010767) q[2];
sx q[2];
rz(-2.0862759) q[2];
rz(-2.1999551) q[3];
sx q[3];
rz(-2.0703273) q[3];
sx q[3];
rz(2.1987703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68614352) q[0];
sx q[0];
rz(-2.2930155) q[0];
sx q[0];
rz(0.69547478) q[0];
rz(-1.6740359) q[1];
sx q[1];
rz(-1.2662042) q[1];
sx q[1];
rz(0.094559018) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3723264) q[0];
sx q[0];
rz(-0.28217286) q[0];
sx q[0];
rz(1.7442987) q[0];
rz(-pi) q[1];
rz(0.90717051) q[2];
sx q[2];
rz(-0.84270515) q[2];
sx q[2];
rz(-1.0502732) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8670204) q[1];
sx q[1];
rz(-2.7944075) q[1];
sx q[1];
rz(-0.95496655) q[1];
rz(-2.2587682) q[3];
sx q[3];
rz(-1.3677639) q[3];
sx q[3];
rz(-0.95850755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5203984) q[2];
sx q[2];
rz(-2.3475671) q[2];
sx q[2];
rz(2.3801079) q[2];
rz(-2.9705808) q[3];
sx q[3];
rz(-1.2820425) q[3];
sx q[3];
rz(0.096693501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.3068202) q[0];
sx q[0];
rz(-2.4691395) q[0];
sx q[0];
rz(-2.6493454) q[0];
rz(-1.3193839) q[1];
sx q[1];
rz(-1.7183813) q[1];
sx q[1];
rz(-0.48670235) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20535417) q[0];
sx q[0];
rz(-1.9689318) q[0];
sx q[0];
rz(-1.4411082) q[0];
rz(-pi) q[1];
rz(0.38551824) q[2];
sx q[2];
rz(-0.80759128) q[2];
sx q[2];
rz(0.68840068) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.69667212) q[1];
sx q[1];
rz(-1.9334403) q[1];
sx q[1];
rz(-1.1128694) q[1];
x q[2];
rz(-0.23426849) q[3];
sx q[3];
rz(-1.0674914) q[3];
sx q[3];
rz(2.4960394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1897366) q[2];
sx q[2];
rz(-2.2978013) q[2];
sx q[2];
rz(2.0437415) q[2];
rz(2.3073933) q[3];
sx q[3];
rz(-2.4662377) q[3];
sx q[3];
rz(1.7043017) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6247691) q[0];
sx q[0];
rz(-2.5747445) q[0];
sx q[0];
rz(-3.095742) q[0];
rz(-0.74371964) q[1];
sx q[1];
rz(-2.7727978) q[1];
sx q[1];
rz(-3.0812982) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9550388) q[0];
sx q[0];
rz(-2.5738724) q[0];
sx q[0];
rz(0.7132775) q[0];
x q[1];
rz(-2.7457775) q[2];
sx q[2];
rz(-0.91378736) q[2];
sx q[2];
rz(2.3558877) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.454596) q[1];
sx q[1];
rz(-2.7424062) q[1];
sx q[1];
rz(2.1737425) q[1];
rz(3.0457959) q[3];
sx q[3];
rz(-2.2471273) q[3];
sx q[3];
rz(-2.6696882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.22659773) q[2];
sx q[2];
rz(-1.3264341) q[2];
sx q[2];
rz(-0.9787406) q[2];
rz(-1.8799051) q[3];
sx q[3];
rz(-2.1224969) q[3];
sx q[3];
rz(-2.252388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9642445) q[0];
sx q[0];
rz(-1.6151936) q[0];
sx q[0];
rz(-0.88441315) q[0];
rz(-2.2824967) q[1];
sx q[1];
rz(-2.0180549) q[1];
sx q[1];
rz(2.9335847) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1387062) q[0];
sx q[0];
rz(-1.0755259) q[0];
sx q[0];
rz(1.1664835) q[0];
x q[1];
rz(-1.154083) q[2];
sx q[2];
rz(-1.0089968) q[2];
sx q[2];
rz(-0.75129189) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0178725) q[1];
sx q[1];
rz(-1.3062638) q[1];
sx q[1];
rz(-0.912028) q[1];
x q[2];
rz(0.0020387928) q[3];
sx q[3];
rz(-1.0674132) q[3];
sx q[3];
rz(0.17755213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7615243) q[2];
sx q[2];
rz(-2.2218349) q[2];
sx q[2];
rz(0.541614) q[2];
rz(-1.8875061) q[3];
sx q[3];
rz(-1.5518291) q[3];
sx q[3];
rz(2.2327173) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33614531) q[0];
sx q[0];
rz(-2.8271524) q[0];
sx q[0];
rz(-1.0910777) q[0];
rz(-2.3167141) q[1];
sx q[1];
rz(-0.42366091) q[1];
sx q[1];
rz(-1.0728015) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0925508) q[0];
sx q[0];
rz(-1.149384) q[0];
sx q[0];
rz(1.3611991) q[0];
rz(-1.5481516) q[2];
sx q[2];
rz(-1.5576819) q[2];
sx q[2];
rz(2.2911521) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8696647) q[1];
sx q[1];
rz(-2.8638693) q[1];
sx q[1];
rz(0.50599338) q[1];
rz(0.24157584) q[3];
sx q[3];
rz(-1.9630057) q[3];
sx q[3];
rz(-1.4763614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.23952809) q[2];
sx q[2];
rz(-2.5967279) q[2];
sx q[2];
rz(5/(6*pi)) q[2];
rz(-0.15051633) q[3];
sx q[3];
rz(-2.033332) q[3];
sx q[3];
rz(0.36293852) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6657669) q[0];
sx q[0];
rz(-3.0463687) q[0];
sx q[0];
rz(1.4775403) q[0];
rz(-2.3382969) q[1];
sx q[1];
rz(-1.7684312) q[1];
sx q[1];
rz(2.1243748) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2870713) q[0];
sx q[0];
rz(-1.6994564) q[0];
sx q[0];
rz(1.3657938) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4847277) q[2];
sx q[2];
rz(-2.629731) q[2];
sx q[2];
rz(-0.32527015) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.89722952) q[1];
sx q[1];
rz(-1.9112226) q[1];
sx q[1];
rz(1.8805518) q[1];
x q[2];
rz(-0.74541645) q[3];
sx q[3];
rz(-2.7099797) q[3];
sx q[3];
rz(3.0109757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9251755) q[2];
sx q[2];
rz(-1.567013) q[2];
sx q[2];
rz(-0.65640059) q[2];
rz(2.9260855) q[3];
sx q[3];
rz(-1.7301205) q[3];
sx q[3];
rz(1.4490674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0553174) q[0];
sx q[0];
rz(-1.7973987) q[0];
sx q[0];
rz(0.84981808) q[0];
rz(-0.96018106) q[1];
sx q[1];
rz(-0.4147059) q[1];
sx q[1];
rz(0.93920952) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16383434) q[0];
sx q[0];
rz(-2.3049666) q[0];
sx q[0];
rz(2.7833582) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1562111) q[2];
sx q[2];
rz(-1.5636946) q[2];
sx q[2];
rz(-2.9455001) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2823742) q[1];
sx q[1];
rz(-1.9422429) q[1];
sx q[1];
rz(-0.66684315) q[1];
x q[2];
rz(-1.0435095) q[3];
sx q[3];
rz(-1.1498972) q[3];
sx q[3];
rz(-2.7038006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5822997) q[2];
sx q[2];
rz(-1.4037932) q[2];
sx q[2];
rz(2.5913008) q[2];
rz(0.127921) q[3];
sx q[3];
rz(-0.13851276) q[3];
sx q[3];
rz(-2.1730455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48548231) q[0];
sx q[0];
rz(-0.67018296) q[0];
sx q[0];
rz(-0.67076587) q[0];
rz(0.35519629) q[1];
sx q[1];
rz(-1.6658446) q[1];
sx q[1];
rz(2.7459941) q[1];
rz(1.4305461) q[2];
sx q[2];
rz(-1.7127019) q[2];
sx q[2];
rz(0.4455926) q[2];
rz(-2.0215423) q[3];
sx q[3];
rz(-1.7786296) q[3];
sx q[3];
rz(2.8854388) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
