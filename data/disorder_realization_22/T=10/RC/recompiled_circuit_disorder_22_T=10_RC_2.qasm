OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7093631) q[0];
sx q[0];
rz(-2.1837283) q[0];
sx q[0];
rz(-0.14444484) q[0];
rz(0.56675178) q[1];
sx q[1];
rz(-0.52539879) q[1];
sx q[1];
rz(0.97775835) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.190783) q[0];
sx q[0];
rz(-1.8291744) q[0];
sx q[0];
rz(1.6590614) q[0];
rz(-2.0830886) q[2];
sx q[2];
rz(-2.5663178) q[2];
sx q[2];
rz(0.011205999) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6482918) q[1];
sx q[1];
rz(-2.9420605) q[1];
sx q[1];
rz(-0.7944016) q[1];
rz(0.17840673) q[3];
sx q[3];
rz(-0.79531407) q[3];
sx q[3];
rz(0.91245302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4686761) q[2];
sx q[2];
rz(-1.1489931) q[2];
sx q[2];
rz(-2.2093175) q[2];
rz(2.9428234) q[3];
sx q[3];
rz(-1.1105744) q[3];
sx q[3];
rz(0.96536243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7070049) q[0];
sx q[0];
rz(-2.2362237) q[0];
sx q[0];
rz(-2.7804651) q[0];
rz(-1.4350285) q[1];
sx q[1];
rz(-1.3577434) q[1];
sx q[1];
rz(-2.3235869) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6011071) q[0];
sx q[0];
rz(-1.2955106) q[0];
sx q[0];
rz(-2.7406373) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8774162) q[2];
sx q[2];
rz(-0.38796705) q[2];
sx q[2];
rz(1.0248794) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.86707838) q[1];
sx q[1];
rz(-2.3550905) q[1];
sx q[1];
rz(2.2134476) q[1];
rz(-1.7198256) q[3];
sx q[3];
rz(-0.36747284) q[3];
sx q[3];
rz(0.056882337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.25257418) q[2];
sx q[2];
rz(-0.44712862) q[2];
sx q[2];
rz(-1.7209631) q[2];
rz(-1.3160926) q[3];
sx q[3];
rz(-2.3836453) q[3];
sx q[3];
rz(-0.38823286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4017568) q[0];
sx q[0];
rz(-0.49210423) q[0];
sx q[0];
rz(0.87093583) q[0];
rz(0.31618205) q[1];
sx q[1];
rz(-0.28156391) q[1];
sx q[1];
rz(-2.8443764) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7446211) q[0];
sx q[0];
rz(-0.40980761) q[0];
sx q[0];
rz(2.9817392) q[0];
rz(-pi) q[1];
rz(-2.259953) q[2];
sx q[2];
rz(-0.92703968) q[2];
sx q[2];
rz(-1.8635441) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.45589046) q[1];
sx q[1];
rz(-2.9487902) q[1];
sx q[1];
rz(1.8222068) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.475297) q[3];
sx q[3];
rz(-0.94216457) q[3];
sx q[3];
rz(0.48376885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7685984) q[2];
sx q[2];
rz(-2.3186389) q[2];
sx q[2];
rz(-2.7139943) q[2];
rz(-1.1887431) q[3];
sx q[3];
rz(-0.62994981) q[3];
sx q[3];
rz(-2.6141613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.4363842) q[0];
sx q[0];
rz(-1.5642865) q[0];
sx q[0];
rz(-0.77392459) q[0];
rz(-2.4286843) q[1];
sx q[1];
rz(-1.0293101) q[1];
sx q[1];
rz(0.68177044) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49127689) q[0];
sx q[0];
rz(-2.804545) q[0];
sx q[0];
rz(0.89462535) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1749928) q[2];
sx q[2];
rz(-1.3912429) q[2];
sx q[2];
rz(2.0008848) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7766429) q[1];
sx q[1];
rz(-0.72307359) q[1];
sx q[1];
rz(-2.2605881) q[1];
rz(-pi) q[2];
x q[2];
rz(0.81687974) q[3];
sx q[3];
rz(-2.6602392) q[3];
sx q[3];
rz(-0.64827418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.008808) q[2];
sx q[2];
rz(-2.4286353) q[2];
sx q[2];
rz(0.95820367) q[2];
rz(-2.0751674) q[3];
sx q[3];
rz(-1.2833779) q[3];
sx q[3];
rz(0.3796033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5565857) q[0];
sx q[0];
rz(-1.3469232) q[0];
sx q[0];
rz(-0.56030309) q[0];
rz(2.141748) q[1];
sx q[1];
rz(-0.20345774) q[1];
sx q[1];
rz(-1.6220185) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9902089) q[0];
sx q[0];
rz(-0.68400331) q[0];
sx q[0];
rz(-0.74599501) q[0];
rz(-pi) q[1];
rz(-2.1278473) q[2];
sx q[2];
rz(-1.0336913) q[2];
sx q[2];
rz(0.86442664) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6474364) q[1];
sx q[1];
rz(-1.7687706) q[1];
sx q[1];
rz(0.76311771) q[1];
rz(-2.1836957) q[3];
sx q[3];
rz(-1.4688204) q[3];
sx q[3];
rz(-1.0979872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6891629) q[2];
sx q[2];
rz(-2.5897557) q[2];
sx q[2];
rz(-2.9186644) q[2];
rz(0.034742268) q[3];
sx q[3];
rz(-1.3870753) q[3];
sx q[3];
rz(0.071578659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8054304) q[0];
sx q[0];
rz(-0.31328377) q[0];
sx q[0];
rz(2.0715332) q[0];
rz(-1.3609715) q[1];
sx q[1];
rz(-0.37934163) q[1];
sx q[1];
rz(-1.8575352) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5491935) q[0];
sx q[0];
rz(-2.1201029) q[0];
sx q[0];
rz(-0.71676371) q[0];
x q[1];
rz(1.9610923) q[2];
sx q[2];
rz(-1.830606) q[2];
sx q[2];
rz(-0.63894546) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9027054) q[1];
sx q[1];
rz(-2.0327912) q[1];
sx q[1];
rz(-1.2483031) q[1];
rz(-pi) q[2];
rz(-2.0606023) q[3];
sx q[3];
rz(-0.59270699) q[3];
sx q[3];
rz(2.5268775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6665035) q[2];
sx q[2];
rz(-1.44375) q[2];
sx q[2];
rz(-0.63759032) q[2];
rz(-2.3049138) q[3];
sx q[3];
rz(-1.1058608) q[3];
sx q[3];
rz(1.3440514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5752983) q[0];
sx q[0];
rz(-0.42995444) q[0];
sx q[0];
rz(-0.56754011) q[0];
rz(0.42770806) q[1];
sx q[1];
rz(-1.6141012) q[1];
sx q[1];
rz(-2.2033851) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8021916) q[0];
sx q[0];
rz(-1.2656478) q[0];
sx q[0];
rz(-1.1423654) q[0];
x q[1];
rz(-2.1392518) q[2];
sx q[2];
rz(-2.217514) q[2];
sx q[2];
rz(1.2127753) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5349622) q[1];
sx q[1];
rz(-2.3772847) q[1];
sx q[1];
rz(2.0984089) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70721831) q[3];
sx q[3];
rz(-1.8054609) q[3];
sx q[3];
rz(-0.16682391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87970916) q[2];
sx q[2];
rz(-1.1738913) q[2];
sx q[2];
rz(-1.3860469) q[2];
rz(-1.8188247) q[3];
sx q[3];
rz(-1.1498007) q[3];
sx q[3];
rz(-3.1125606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8742074) q[0];
sx q[0];
rz(-2.8081942) q[0];
sx q[0];
rz(1.4338795) q[0];
rz(1.2738312) q[1];
sx q[1];
rz(-1.1359943) q[1];
sx q[1];
rz(-2.3103255) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2682876) q[0];
sx q[0];
rz(-1.9718861) q[0];
sx q[0];
rz(1.7307161) q[0];
x q[1];
rz(-1.7449042) q[2];
sx q[2];
rz(-0.6422407) q[2];
sx q[2];
rz(1.0345392) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9030994) q[1];
sx q[1];
rz(-1.0817173) q[1];
sx q[1];
rz(-0.78519435) q[1];
x q[2];
rz(1.3475111) q[3];
sx q[3];
rz(-1.7956453) q[3];
sx q[3];
rz(-2.4527578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5635809) q[2];
sx q[2];
rz(-1.3728023) q[2];
sx q[2];
rz(-0.9643628) q[2];
rz(-1.1635121) q[3];
sx q[3];
rz(-0.61146277) q[3];
sx q[3];
rz(-0.65892974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7678541) q[0];
sx q[0];
rz(-2.4066194) q[0];
sx q[0];
rz(2.1642165) q[0];
rz(1.7550229) q[1];
sx q[1];
rz(-1.8354548) q[1];
sx q[1];
rz(-2.0358553) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10406636) q[0];
sx q[0];
rz(-0.46188799) q[0];
sx q[0];
rz(0.26432963) q[0];
x q[1];
rz(0.10496225) q[2];
sx q[2];
rz(-2.2347921) q[2];
sx q[2];
rz(2.2763989) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2178206) q[1];
sx q[1];
rz(-2.8433228) q[1];
sx q[1];
rz(-2.24733) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4615199) q[3];
sx q[3];
rz(-2.701093) q[3];
sx q[3];
rz(-2.3624453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5332807) q[2];
sx q[2];
rz(-2.7567342) q[2];
sx q[2];
rz(-2.9917955) q[2];
rz(1.3730565) q[3];
sx q[3];
rz(-1.7356197) q[3];
sx q[3];
rz(1.3214553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49366429) q[0];
sx q[0];
rz(-2.6239008) q[0];
sx q[0];
rz(0.1272442) q[0];
rz(1.6607025) q[1];
sx q[1];
rz(-2.7791185) q[1];
sx q[1];
rz(0.1677992) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0566933) q[0];
sx q[0];
rz(-3.1057682) q[0];
sx q[0];
rz(-2.5762659) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5334353) q[2];
sx q[2];
rz(-1.4434575) q[2];
sx q[2];
rz(-1.5978158) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2192981) q[1];
sx q[1];
rz(-2.1592327) q[1];
sx q[1];
rz(1.5625619) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76370244) q[3];
sx q[3];
rz(-2.9716431) q[3];
sx q[3];
rz(-0.22596879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.026713513) q[2];
sx q[2];
rz(-2.202704) q[2];
sx q[2];
rz(-0.76114571) q[2];
rz(0.090027697) q[3];
sx q[3];
rz(-2.138425) q[3];
sx q[3];
rz(0.95054039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54820838) q[0];
sx q[0];
rz(-1.1593288) q[0];
sx q[0];
rz(2.819084) q[0];
rz(-2.7535915) q[1];
sx q[1];
rz(-1.7419659) q[1];
sx q[1];
rz(2.3566125) q[1];
rz(-0.81007304) q[2];
sx q[2];
rz(-2.0740866) q[2];
sx q[2];
rz(-1.496051) q[2];
rz(2.1098059) q[3];
sx q[3];
rz(-1.7179334) q[3];
sx q[3];
rz(-2.4739305) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
