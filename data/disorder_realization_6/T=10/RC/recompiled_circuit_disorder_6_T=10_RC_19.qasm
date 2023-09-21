OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6157827) q[0];
sx q[0];
rz(-1.4178185) q[0];
sx q[0];
rz(0.56086993) q[0];
rz(1.1129192) q[1];
sx q[1];
rz(-1.7634044) q[1];
sx q[1];
rz(1.2150432) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96072223) q[0];
sx q[0];
rz(-1.3836432) q[0];
sx q[0];
rz(0.10589177) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66081337) q[2];
sx q[2];
rz(-1.8006969) q[2];
sx q[2];
rz(1.4290609) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3455968) q[1];
sx q[1];
rz(-1.0618292) q[1];
sx q[1];
rz(0.3791581) q[1];
rz(-pi) q[2];
rz(1.8362942) q[3];
sx q[3];
rz(-1.7146829) q[3];
sx q[3];
rz(-0.063751566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3540196) q[2];
sx q[2];
rz(-2.188787) q[2];
sx q[2];
rz(-0.18307486) q[2];
rz(-0.37781528) q[3];
sx q[3];
rz(-1.0487882) q[3];
sx q[3];
rz(-2.8474076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8437682) q[0];
sx q[0];
rz(-2.4968708) q[0];
sx q[0];
rz(0.077117292) q[0];
rz(2.8027957) q[1];
sx q[1];
rz(-2.0270551) q[1];
sx q[1];
rz(-1.6024626) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3829271) q[0];
sx q[0];
rz(-0.59249741) q[0];
sx q[0];
rz(-1.7090319) q[0];
x q[1];
rz(1.7878754) q[2];
sx q[2];
rz(-2.2976544) q[2];
sx q[2];
rz(0.94490563) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7669249) q[1];
sx q[1];
rz(-2.6267849) q[1];
sx q[1];
rz(1.2109846) q[1];
rz(-pi) q[2];
rz(0.48861309) q[3];
sx q[3];
rz(-2.6693137) q[3];
sx q[3];
rz(0.74913914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2960647) q[2];
sx q[2];
rz(-1.277593) q[2];
sx q[2];
rz(2.4831333) q[2];
rz(-0.15130875) q[3];
sx q[3];
rz(-2.1189809) q[3];
sx q[3];
rz(-0.69491274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.113134) q[0];
sx q[0];
rz(-2.3582393) q[0];
sx q[0];
rz(-0.43310305) q[0];
rz(-1.1921047) q[1];
sx q[1];
rz(-1.9299709) q[1];
sx q[1];
rz(-0.55535299) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1028324) q[0];
sx q[0];
rz(-2.0534678) q[0];
sx q[0];
rz(-2.32248) q[0];
rz(1.0054587) q[2];
sx q[2];
rz(-1.8381422) q[2];
sx q[2];
rz(2.8440059) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.66545031) q[1];
sx q[1];
rz(-0.94201554) q[1];
sx q[1];
rz(2.9412342) q[1];
rz(-1.2461353) q[3];
sx q[3];
rz(-2.5723296) q[3];
sx q[3];
rz(0.07490052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.73734036) q[2];
sx q[2];
rz(-2.3534687) q[2];
sx q[2];
rz(1.2505442) q[2];
rz(-0.2441497) q[3];
sx q[3];
rz(-1.8593676) q[3];
sx q[3];
rz(-1.6916493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26043949) q[0];
sx q[0];
rz(-0.45757159) q[0];
sx q[0];
rz(0.81480169) q[0];
rz(1.762215) q[1];
sx q[1];
rz(-0.35019362) q[1];
sx q[1];
rz(0.25517685) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4677306) q[0];
sx q[0];
rz(-2.6669589) q[0];
sx q[0];
rz(-0.69068308) q[0];
x q[1];
rz(2.5023979) q[2];
sx q[2];
rz(-1.3022458) q[2];
sx q[2];
rz(-0.07721363) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0604917) q[1];
sx q[1];
rz(-1.2730518) q[1];
sx q[1];
rz(-0.22025073) q[1];
rz(-pi) q[2];
rz(1.4820443) q[3];
sx q[3];
rz(-2.6188861) q[3];
sx q[3];
rz(-2.6356217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.8884376) q[2];
sx q[2];
rz(-1.5909114) q[2];
sx q[2];
rz(2.9684084) q[2];
rz(-2.611768) q[3];
sx q[3];
rz(-0.14557043) q[3];
sx q[3];
rz(3.0392652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2816876) q[0];
sx q[0];
rz(-1.4929993) q[0];
sx q[0];
rz(1.3758855) q[0];
rz(1.2777404) q[1];
sx q[1];
rz(-2.3294096) q[1];
sx q[1];
rz(-3.0854991) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027095196) q[0];
sx q[0];
rz(-1.9341015) q[0];
sx q[0];
rz(-2.5278805) q[0];
rz(-pi) q[1];
rz(-1.1115083) q[2];
sx q[2];
rz(-1.2570724) q[2];
sx q[2];
rz(0.26278652) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5406815) q[1];
sx q[1];
rz(-1.8735421) q[1];
sx q[1];
rz(-1.5860228) q[1];
rz(-pi) q[2];
rz(-1.7124743) q[3];
sx q[3];
rz(-0.56483993) q[3];
sx q[3];
rz(-0.39839881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1405979) q[2];
sx q[2];
rz(-2.2237015) q[2];
sx q[2];
rz(-0.43219217) q[2];
rz(2.2473992) q[3];
sx q[3];
rz(-1.0995355) q[3];
sx q[3];
rz(1.6754707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11319259) q[0];
sx q[0];
rz(-2.2560461) q[0];
sx q[0];
rz(2.4940441) q[0];
rz(1.8796857) q[1];
sx q[1];
rz(-1.6779265) q[1];
sx q[1];
rz(-2.1870959) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0237085) q[0];
sx q[0];
rz(-0.55340289) q[0];
sx q[0];
rz(1.2721567) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1068222) q[2];
sx q[2];
rz(-2.2010942) q[2];
sx q[2];
rz(-2.6882753) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.33659014) q[1];
sx q[1];
rz(-1.0953566) q[1];
sx q[1];
rz(-2.5617983) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0733842) q[3];
sx q[3];
rz(-0.32940255) q[3];
sx q[3];
rz(-2.9941032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.548617) q[2];
sx q[2];
rz(-1.9081215) q[2];
sx q[2];
rz(-2.0992289) q[2];
rz(2.7029165) q[3];
sx q[3];
rz(-2.091566) q[3];
sx q[3];
rz(1.3180805) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2838659) q[0];
sx q[0];
rz(-0.23290578) q[0];
sx q[0];
rz(2.3983811) q[0];
rz(-1.5076393) q[1];
sx q[1];
rz(-2.4217024) q[1];
sx q[1];
rz(-2.5315703) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7689432) q[0];
sx q[0];
rz(-0.56108755) q[0];
sx q[0];
rz(2.6555496) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0257752) q[2];
sx q[2];
rz(-2.2216185) q[2];
sx q[2];
rz(1.7388369) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4808018) q[1];
sx q[1];
rz(-1.5958438) q[1];
sx q[1];
rz(-2.2381496) q[1];
rz(-pi) q[2];
rz(-1.0878629) q[3];
sx q[3];
rz(-2.7186839) q[3];
sx q[3];
rz(1.9382697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8053749) q[2];
sx q[2];
rz(-1.6992133) q[2];
sx q[2];
rz(0.91840333) q[2];
rz(1.5911128) q[3];
sx q[3];
rz(-2.1912626) q[3];
sx q[3];
rz(-0.38890719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36088762) q[0];
sx q[0];
rz(-0.66910678) q[0];
sx q[0];
rz(1.6280744) q[0];
rz(2.6121415) q[1];
sx q[1];
rz(-1.0667195) q[1];
sx q[1];
rz(-2.4050074) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1153699) q[0];
sx q[0];
rz(-1.5407469) q[0];
sx q[0];
rz(-3.1307334) q[0];
rz(0.10891624) q[2];
sx q[2];
rz(-1.8913336) q[2];
sx q[2];
rz(-0.045217302) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8383933) q[1];
sx q[1];
rz(-2.0337078) q[1];
sx q[1];
rz(-1.4569015) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1804817) q[3];
sx q[3];
rz(-2.612252) q[3];
sx q[3];
rz(1.6314268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0344051) q[2];
sx q[2];
rz(-1.2074869) q[2];
sx q[2];
rz(0.68391189) q[2];
rz(-1.9125787) q[3];
sx q[3];
rz(-1.3701655) q[3];
sx q[3];
rz(-1.3945403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1972315) q[0];
sx q[0];
rz(-1.5427538) q[0];
sx q[0];
rz(2.9558682) q[0];
rz(-2.1445403) q[1];
sx q[1];
rz(-1.2652218) q[1];
sx q[1];
rz(2.396778) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54430994) q[0];
sx q[0];
rz(-2.3047857) q[0];
sx q[0];
rz(0.44041667) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3486175) q[2];
sx q[2];
rz(-1.8599531) q[2];
sx q[2];
rz(2.5536429) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3265171) q[1];
sx q[1];
rz(-2.7217367) q[1];
sx q[1];
rz(2.144787) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7471077) q[3];
sx q[3];
rz(-2.3186765) q[3];
sx q[3];
rz(-2.8701973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0361438) q[2];
sx q[2];
rz(-0.75573409) q[2];
sx q[2];
rz(-1.9469117) q[2];
rz(-2.1448994) q[3];
sx q[3];
rz(-1.9255305) q[3];
sx q[3];
rz(-0.99635807) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74334082) q[0];
sx q[0];
rz(-2.353459) q[0];
sx q[0];
rz(-0.40400305) q[0];
rz(-0.031127302) q[1];
sx q[1];
rz(-1.6571836) q[1];
sx q[1];
rz(1.1709447) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4963213) q[0];
sx q[0];
rz(-1.5924581) q[0];
sx q[0];
rz(2.7146656) q[0];
x q[1];
rz(-1.3002214) q[2];
sx q[2];
rz(-2.3279394) q[2];
sx q[2];
rz(2.4547581) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.4198235) q[1];
sx q[1];
rz(-3.1173919) q[1];
sx q[1];
rz(-1.9355378) q[1];
x q[2];
rz(1.8627432) q[3];
sx q[3];
rz(-1.1598831) q[3];
sx q[3];
rz(1.7850072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4460454) q[2];
sx q[2];
rz(-1.7798767) q[2];
sx q[2];
rz(-0.5919624) q[2];
rz(-2.5752318) q[3];
sx q[3];
rz(-0.16470328) q[3];
sx q[3];
rz(-1.6177572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82407172) q[0];
sx q[0];
rz(-0.98012797) q[0];
sx q[0];
rz(-1.160887) q[0];
rz(-0.099427632) q[1];
sx q[1];
rz(-1.2482523) q[1];
sx q[1];
rz(-2.0773239) q[1];
rz(2.4026985) q[2];
sx q[2];
rz(-2.0901491) q[2];
sx q[2];
rz(0.73170589) q[2];
rz(-0.016146544) q[3];
sx q[3];
rz(-1.2373677) q[3];
sx q[3];
rz(2.2131372) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];