OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.82436615) q[0];
sx q[0];
rz(-1.1146201) q[0];
sx q[0];
rz(-0.00014076509) q[0];
rz(1.3340985) q[1];
sx q[1];
rz(4.1058022) q[1];
sx q[1];
rz(10.618187) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1228468) q[0];
sx q[0];
rz(-1.7483286) q[0];
sx q[0];
rz(-1.2794504) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.869007) q[2];
sx q[2];
rz(-1.0422921) q[2];
sx q[2];
rz(-0.82982153) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.57230091) q[1];
sx q[1];
rz(-0.83218677) q[1];
sx q[1];
rz(-2.4898847) q[1];
rz(-pi) q[2];
rz(-1.1080997) q[3];
sx q[3];
rz(-0.24216147) q[3];
sx q[3];
rz(-2.7778181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.45941916) q[2];
sx q[2];
rz(-0.023962263) q[2];
sx q[2];
rz(-1.2288644) q[2];
rz(-1.7284283) q[3];
sx q[3];
rz(-2.0404405) q[3];
sx q[3];
rz(-1.4878954) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6035778) q[0];
sx q[0];
rz(-1.5025654) q[0];
sx q[0];
rz(2.1287825) q[0];
rz(3.1139328) q[1];
sx q[1];
rz(-2.467997) q[1];
sx q[1];
rz(2.0181296) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2259953) q[0];
sx q[0];
rz(-1.628327) q[0];
sx q[0];
rz(2.0321839) q[0];
x q[1];
rz(0.77983071) q[2];
sx q[2];
rz(-1.6155598) q[2];
sx q[2];
rz(-1.1080527) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0151129) q[1];
sx q[1];
rz(-2.0868595) q[1];
sx q[1];
rz(-2.5491787) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4853737) q[3];
sx q[3];
rz(-1.1074293) q[3];
sx q[3];
rz(3.0955293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.79364395) q[2];
sx q[2];
rz(-2.0517893) q[2];
sx q[2];
rz(2.222555) q[2];
rz(0.67409003) q[3];
sx q[3];
rz(-0.6522817) q[3];
sx q[3];
rz(-1.526171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27750257) q[0];
sx q[0];
rz(-0.16177495) q[0];
sx q[0];
rz(1.8664237) q[0];
rz(-2.4480942) q[1];
sx q[1];
rz(-1.8854515) q[1];
sx q[1];
rz(-1.1330053) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64273524) q[0];
sx q[0];
rz(-0.14232902) q[0];
sx q[0];
rz(-1.5940773) q[0];
rz(-pi) q[1];
rz(0.79046952) q[2];
sx q[2];
rz(-1.9445473) q[2];
sx q[2];
rz(-0.30465301) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9085711) q[1];
sx q[1];
rz(-2.1070707) q[1];
sx q[1];
rz(0.33645333) q[1];
rz(-0.95136178) q[3];
sx q[3];
rz(-2.1054483) q[3];
sx q[3];
rz(1.9922647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8901849) q[2];
sx q[2];
rz(-0.79139411) q[2];
sx q[2];
rz(1.8481002) q[2];
rz(-3.1022762) q[3];
sx q[3];
rz(-1.2189564) q[3];
sx q[3];
rz(1.2600651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2599729) q[0];
sx q[0];
rz(-3.0631174) q[0];
sx q[0];
rz(1.1608634) q[0];
rz(-0.89598957) q[1];
sx q[1];
rz(-1.4410102) q[1];
sx q[1];
rz(3.0060351) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7318152) q[0];
sx q[0];
rz(-1.2588132) q[0];
sx q[0];
rz(0.72809763) q[0];
x q[1];
rz(1.3877669) q[2];
sx q[2];
rz(-0.39003885) q[2];
sx q[2];
rz(-2.4412465) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.902066) q[1];
sx q[1];
rz(-2.8956928) q[1];
sx q[1];
rz(-2.2794101) q[1];
rz(-1.0566063) q[3];
sx q[3];
rz(-2.8286472) q[3];
sx q[3];
rz(1.3514951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23665145) q[2];
sx q[2];
rz(-0.94649482) q[2];
sx q[2];
rz(2.2616852) q[2];
rz(3.0974292) q[3];
sx q[3];
rz(-1.6396089) q[3];
sx q[3];
rz(-2.8529609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0376461) q[0];
sx q[0];
rz(-0.3750616) q[0];
sx q[0];
rz(-2.1283545) q[0];
rz(3.0918616) q[1];
sx q[1];
rz(-0.91369349) q[1];
sx q[1];
rz(2.0577046) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.097682) q[0];
sx q[0];
rz(-1.327276) q[0];
sx q[0];
rz(-2.9527412) q[0];
rz(-0.25123698) q[2];
sx q[2];
rz(-1.8357956) q[2];
sx q[2];
rz(2.3064409) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.14099546) q[1];
sx q[1];
rz(-0.53783572) q[1];
sx q[1];
rz(-1.7829893) q[1];
rz(0.089133457) q[3];
sx q[3];
rz(-0.51557487) q[3];
sx q[3];
rz(-2.6775132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.23285398) q[2];
sx q[2];
rz(-2.8149657) q[2];
sx q[2];
rz(0.24442913) q[2];
rz(-2.7092253) q[3];
sx q[3];
rz(-1.3997388) q[3];
sx q[3];
rz(-0.50306815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8571092) q[0];
sx q[0];
rz(-1.720022) q[0];
sx q[0];
rz(-3.0474512) q[0];
rz(-0.17177467) q[1];
sx q[1];
rz(-2.005902) q[1];
sx q[1];
rz(-2.24618) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58127922) q[0];
sx q[0];
rz(-0.34438294) q[0];
sx q[0];
rz(-3.0292105) q[0];
rz(-pi) q[1];
rz(2.5080639) q[2];
sx q[2];
rz(-1.3760566) q[2];
sx q[2];
rz(-0.63846904) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.15570607) q[1];
sx q[1];
rz(-2.967318) q[1];
sx q[1];
rz(1.8272912) q[1];
rz(-pi) q[2];
rz(1.5395245) q[3];
sx q[3];
rz(-0.074723738) q[3];
sx q[3];
rz(-0.63262313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0078997) q[2];
sx q[2];
rz(-0.40863016) q[2];
sx q[2];
rz(-0.80319476) q[2];
rz(1.1903654) q[3];
sx q[3];
rz(-1.9093711) q[3];
sx q[3];
rz(-2.7289594) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0728834) q[0];
sx q[0];
rz(-2.9769653) q[0];
sx q[0];
rz(2.6224526) q[0];
rz(-0.58147645) q[1];
sx q[1];
rz(-1.1053718) q[1];
sx q[1];
rz(1.2566459) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089766895) q[0];
sx q[0];
rz(-1.8221812) q[0];
sx q[0];
rz(0.041640394) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2101375) q[2];
sx q[2];
rz(-1.1147095) q[2];
sx q[2];
rz(1.2880585) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.41374846) q[1];
sx q[1];
rz(-1.2894948) q[1];
sx q[1];
rz(-1.6378228) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.890896) q[3];
sx q[3];
rz(-2.2543636) q[3];
sx q[3];
rz(-0.86910955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1533623) q[2];
sx q[2];
rz(-1.0299269) q[2];
sx q[2];
rz(1.777565) q[2];
rz(-2.2310232) q[3];
sx q[3];
rz(-1.1547337) q[3];
sx q[3];
rz(1.6114657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0664739) q[0];
sx q[0];
rz(-2.5771038) q[0];
sx q[0];
rz(0.30817729) q[0];
rz(0.072487436) q[1];
sx q[1];
rz(-1.0132353) q[1];
sx q[1];
rz(2.7546308) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.455084) q[0];
sx q[0];
rz(-0.63502705) q[0];
sx q[0];
rz(-1.9576661) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96078028) q[2];
sx q[2];
rz(-0.79723251) q[2];
sx q[2];
rz(-1.0459935) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9427529) q[1];
sx q[1];
rz(-1.6457335) q[1];
sx q[1];
rz(-2.4692175) q[1];
rz(-pi) q[2];
rz(-2.7510795) q[3];
sx q[3];
rz(-0.53005855) q[3];
sx q[3];
rz(0.9583677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5179634) q[2];
sx q[2];
rz(-1.364418) q[2];
sx q[2];
rz(2.7015838) q[2];
rz(0.7157588) q[3];
sx q[3];
rz(-1.4322759) q[3];
sx q[3];
rz(1.0796775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14116645) q[0];
sx q[0];
rz(-0.74581242) q[0];
sx q[0];
rz(-1.0986885) q[0];
rz(-2.4138342) q[1];
sx q[1];
rz(-0.37574238) q[1];
sx q[1];
rz(-0.049302014) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44137529) q[0];
sx q[0];
rz(-0.8964296) q[0];
sx q[0];
rz(-0.81654878) q[0];
rz(-pi) q[1];
rz(0.88843139) q[2];
sx q[2];
rz(-1.1738136) q[2];
sx q[2];
rz(1.9156485) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8763435) q[1];
sx q[1];
rz(-2.0298404) q[1];
sx q[1];
rz(-2.6190119) q[1];
rz(-pi) q[2];
rz(-0.74069549) q[3];
sx q[3];
rz(-1.1974088) q[3];
sx q[3];
rz(1.545056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0631642) q[2];
sx q[2];
rz(-0.59946632) q[2];
sx q[2];
rz(-2.4196529) q[2];
rz(2.1980964) q[3];
sx q[3];
rz(-2.3908581) q[3];
sx q[3];
rz(0.25434428) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14770517) q[0];
sx q[0];
rz(-1.9839956) q[0];
sx q[0];
rz(2.06185) q[0];
rz(-2.0823157) q[1];
sx q[1];
rz(-0.22288999) q[1];
sx q[1];
rz(1.7396897) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3831543) q[0];
sx q[0];
rz(-0.24807319) q[0];
sx q[0];
rz(-2.1092578) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5435113) q[2];
sx q[2];
rz(-1.7575043) q[2];
sx q[2];
rz(0.4208496) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.68114963) q[1];
sx q[1];
rz(-2.3773758) q[1];
sx q[1];
rz(2.0185673) q[1];
rz(-pi) q[2];
rz(2.510342) q[3];
sx q[3];
rz(-1.7389986) q[3];
sx q[3];
rz(1.2390618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6587276) q[2];
sx q[2];
rz(-1.3663224) q[2];
sx q[2];
rz(-1.520291) q[2];
rz(2.5907717) q[3];
sx q[3];
rz(-2.3362624) q[3];
sx q[3];
rz(-0.6974535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.993492) q[0];
sx q[0];
rz(-1.3052595) q[0];
sx q[0];
rz(-1.530151) q[0];
rz(0.91611721) q[1];
sx q[1];
rz(-0.59090186) q[1];
sx q[1];
rz(-0.59060243) q[1];
rz(2.312071) q[2];
sx q[2];
rz(-0.98486949) q[2];
sx q[2];
rz(-2.9688901) q[2];
rz(0.017833088) q[3];
sx q[3];
rz(-1.0537536) q[3];
sx q[3];
rz(-0.23585933) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];