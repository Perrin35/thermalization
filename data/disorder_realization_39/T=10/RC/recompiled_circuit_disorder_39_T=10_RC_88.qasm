OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0948148) q[0];
sx q[0];
rz(4.2098213) q[0];
sx q[0];
rz(9.8888483) q[0];
rz(1.9595454) q[1];
sx q[1];
rz(-0.067117604) q[1];
sx q[1];
rz(1.2844515) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5285437) q[0];
sx q[0];
rz(-1.9595946) q[0];
sx q[0];
rz(2.8253428) q[0];
rz(-pi) q[1];
rz(0.061878248) q[2];
sx q[2];
rz(-0.49052325) q[2];
sx q[2];
rz(2.2762736) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1661108) q[1];
sx q[1];
rz(-1.2088641) q[1];
sx q[1];
rz(-2.6544177) q[1];
rz(2.9908882) q[3];
sx q[3];
rz(-1.5871443) q[3];
sx q[3];
rz(-2.5022262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.39711943) q[2];
sx q[2];
rz(-2.2019272) q[2];
sx q[2];
rz(-0.31952566) q[2];
rz(-0.5784936) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(2.4676676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7330866) q[0];
sx q[0];
rz(-1.5717614) q[0];
sx q[0];
rz(2.615036) q[0];
rz(-2.5780442) q[1];
sx q[1];
rz(-1.6105517) q[1];
sx q[1];
rz(2.3449576) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5027673) q[0];
sx q[0];
rz(-0.728038) q[0];
sx q[0];
rz(2.0730221) q[0];
x q[1];
rz(2.8393306) q[2];
sx q[2];
rz(-0.59603359) q[2];
sx q[2];
rz(-0.1421393) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7318774) q[1];
sx q[1];
rz(-2.5341923) q[1];
sx q[1];
rz(-3.1190447) q[1];
rz(-2.6437831) q[3];
sx q[3];
rz(-1.0199162) q[3];
sx q[3];
rz(2.2450972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.18156302) q[2];
sx q[2];
rz(-1.3788362) q[2];
sx q[2];
rz(-0.78655085) q[2];
rz(-0.49318796) q[3];
sx q[3];
rz(-1.1816053) q[3];
sx q[3];
rz(-2.6942159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86984533) q[0];
sx q[0];
rz(-0.87690502) q[0];
sx q[0];
rz(1.440381) q[0];
rz(0.72021833) q[1];
sx q[1];
rz(-0.59599382) q[1];
sx q[1];
rz(-2.4386141) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36613208) q[0];
sx q[0];
rz(-1.6876939) q[0];
sx q[0];
rz(1.2890105) q[0];
rz(-pi) q[1];
rz(1.9269283) q[2];
sx q[2];
rz(-1.0991569) q[2];
sx q[2];
rz(2.4871662) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8382593) q[1];
sx q[1];
rz(-1.6175744) q[1];
sx q[1];
rz(-0.10951885) q[1];
rz(-pi) q[2];
rz(-2.968077) q[3];
sx q[3];
rz(-0.95458889) q[3];
sx q[3];
rz(-0.71615744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.26677033) q[2];
sx q[2];
rz(-1.6569933) q[2];
sx q[2];
rz(-0.91153574) q[2];
rz(-0.91397816) q[3];
sx q[3];
rz(-1.4303047) q[3];
sx q[3];
rz(2.3142557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5754958) q[0];
sx q[0];
rz(-0.63593447) q[0];
sx q[0];
rz(-2.3655868) q[0];
rz(1.2674747) q[1];
sx q[1];
rz(-1.1088561) q[1];
sx q[1];
rz(-0.56328303) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3455428) q[0];
sx q[0];
rz(-0.69957083) q[0];
sx q[0];
rz(-1.1388586) q[0];
rz(2.7411555) q[2];
sx q[2];
rz(-2.4106328) q[2];
sx q[2];
rz(1.8797344) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5717585) q[1];
sx q[1];
rz(-2.8697439) q[1];
sx q[1];
rz(0.024072577) q[1];
rz(-pi) q[2];
rz(-0.33470811) q[3];
sx q[3];
rz(-2.3957806) q[3];
sx q[3];
rz(-2.0832182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6528066) q[2];
sx q[2];
rz(-2.2061009) q[2];
sx q[2];
rz(-0.164786) q[2];
rz(2.9131043) q[3];
sx q[3];
rz(-0.27992862) q[3];
sx q[3];
rz(-2.1951108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27424681) q[0];
sx q[0];
rz(-2.2773401) q[0];
sx q[0];
rz(-0.85246032) q[0];
rz(-2.7903941) q[1];
sx q[1];
rz(-2.5636667) q[1];
sx q[1];
rz(1.9794827) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6974555) q[0];
sx q[0];
rz(-1.9636969) q[0];
sx q[0];
rz(0.59209728) q[0];
x q[1];
rz(-2.2010872) q[2];
sx q[2];
rz(-2.5224707) q[2];
sx q[2];
rz(-1.8096015) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3106766) q[1];
sx q[1];
rz(-1.984593) q[1];
sx q[1];
rz(-0.78891854) q[1];
x q[2];
rz(-0.72158738) q[3];
sx q[3];
rz(-1.5820832) q[3];
sx q[3];
rz(-1.23502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44624415) q[2];
sx q[2];
rz(-0.39351311) q[2];
sx q[2];
rz(1.2472786) q[2];
rz(3.128483) q[3];
sx q[3];
rz(-1.7316827) q[3];
sx q[3];
rz(0.22918992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.2728249) q[0];
sx q[0];
rz(-1.8936963) q[0];
sx q[0];
rz(2.390958) q[0];
rz(-2.0320832) q[1];
sx q[1];
rz(-0.37056071) q[1];
sx q[1];
rz(-1.8355339) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.437498) q[0];
sx q[0];
rz(-1.1878345) q[0];
sx q[0];
rz(-1.0136481) q[0];
rz(1.7734217) q[2];
sx q[2];
rz(-2.0261923) q[2];
sx q[2];
rz(0.95603285) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8522779) q[1];
sx q[1];
rz(-1.5783974) q[1];
sx q[1];
rz(0.72035933) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4218876) q[3];
sx q[3];
rz(-2.1544666) q[3];
sx q[3];
rz(-0.35051171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3874454) q[2];
sx q[2];
rz(-1.0070609) q[2];
sx q[2];
rz(-0.60069096) q[2];
rz(1.0026275) q[3];
sx q[3];
rz(-2.6032344) q[3];
sx q[3];
rz(0.35282648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3951185) q[0];
sx q[0];
rz(-1.4557319) q[0];
sx q[0];
rz(2.0948998) q[0];
rz(1.5294317) q[1];
sx q[1];
rz(-1.4271095) q[1];
sx q[1];
rz(-0.41710645) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0171623) q[0];
sx q[0];
rz(-2.9201047) q[0];
sx q[0];
rz(-1.682196) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.276937) q[2];
sx q[2];
rz(-1.8655348) q[2];
sx q[2];
rz(-0.35640946) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5187832) q[1];
sx q[1];
rz(-0.86474027) q[1];
sx q[1];
rz(2.519636) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.65317513) q[3];
sx q[3];
rz(-1.8110868) q[3];
sx q[3];
rz(-1.0550635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1356915) q[2];
sx q[2];
rz(-0.4726755) q[2];
sx q[2];
rz(1.7741514) q[2];
rz(-2.588429) q[3];
sx q[3];
rz(-1.7382712) q[3];
sx q[3];
rz(0.52545351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.816514) q[0];
sx q[0];
rz(-1.8186318) q[0];
sx q[0];
rz(-1.1897855) q[0];
rz(1.7143543) q[1];
sx q[1];
rz(-0.27806565) q[1];
sx q[1];
rz(0.11238012) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1988797) q[0];
sx q[0];
rz(-2.6330224) q[0];
sx q[0];
rz(-1.2099427) q[0];
rz(-pi) q[1];
rz(-1.4833647) q[2];
sx q[2];
rz(-1.9033252) q[2];
sx q[2];
rz(-2.33193) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5210133) q[1];
sx q[1];
rz(-0.81110209) q[1];
sx q[1];
rz(1.0902507) q[1];
rz(-pi) q[2];
rz(1.8295248) q[3];
sx q[3];
rz(-1.4530164) q[3];
sx q[3];
rz(-0.044101276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0671976) q[2];
sx q[2];
rz(-1.9125568) q[2];
sx q[2];
rz(-1.3486264) q[2];
rz(-1.9366692) q[3];
sx q[3];
rz(-1.7347074) q[3];
sx q[3];
rz(-2.9437734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39712054) q[0];
sx q[0];
rz(-1.4646126) q[0];
sx q[0];
rz(-0.034974139) q[0];
rz(-2.2947625) q[1];
sx q[1];
rz(-1.7085107) q[1];
sx q[1];
rz(0.91167489) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9104886) q[0];
sx q[0];
rz(-1.7945053) q[0];
sx q[0];
rz(1.847812) q[0];
rz(-2.9334925) q[2];
sx q[2];
rz(-0.55317438) q[2];
sx q[2];
rz(2.2262239) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4166491) q[1];
sx q[1];
rz(-2.0282201) q[1];
sx q[1];
rz(-1.7979421) q[1];
rz(-1.0802286) q[3];
sx q[3];
rz(-2.053223) q[3];
sx q[3];
rz(2.3232943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.70790616) q[2];
sx q[2];
rz(-0.5138548) q[2];
sx q[2];
rz(-2.8923477) q[2];
rz(-2.3748659) q[3];
sx q[3];
rz(-1.8079575) q[3];
sx q[3];
rz(0.39961091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3578167) q[0];
sx q[0];
rz(-2.0932842) q[0];
sx q[0];
rz(-2.7695079) q[0];
rz(2.5601939) q[1];
sx q[1];
rz(-0.77722469) q[1];
sx q[1];
rz(1.7262329) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5370731) q[0];
sx q[0];
rz(-1.8329289) q[0];
sx q[0];
rz(-1.7781236) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62016983) q[2];
sx q[2];
rz(-1.1405186) q[2];
sx q[2];
rz(2.534453) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9006173) q[1];
sx q[1];
rz(-1.6458578) q[1];
sx q[1];
rz(1.5013298) q[1];
x q[2];
rz(-2.1595702) q[3];
sx q[3];
rz(-0.62871274) q[3];
sx q[3];
rz(-0.94129896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8156478) q[2];
sx q[2];
rz(-0.14020136) q[2];
sx q[2];
rz(2.9369205) q[2];
rz(-1.7278016) q[3];
sx q[3];
rz(-1.1433583) q[3];
sx q[3];
rz(1.0958825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6407912) q[0];
sx q[0];
rz(-1.4688107) q[0];
sx q[0];
rz(-2.2039913) q[0];
rz(-1.5564556) q[1];
sx q[1];
rz(-1.7592447) q[1];
sx q[1];
rz(-1.6978282) q[1];
rz(1.8150868) q[2];
sx q[2];
rz(-0.46637022) q[2];
sx q[2];
rz(-2.7844219) q[2];
rz(-3.0588991) q[3];
sx q[3];
rz(-2.1413998) q[3];
sx q[3];
rz(-1.0257046) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
