OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2919579) q[0];
sx q[0];
rz(-2.7014974) q[0];
sx q[0];
rz(3.0043998) q[0];
rz(1.4057012) q[1];
sx q[1];
rz(-1.7383716) q[1];
sx q[1];
rz(0.52991968) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0041381) q[0];
sx q[0];
rz(-1.19085) q[0];
sx q[0];
rz(-0.11560346) q[0];
x q[1];
rz(-1.1692739) q[2];
sx q[2];
rz(-1.2520773) q[2];
sx q[2];
rz(-0.67414325) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.73218988) q[1];
sx q[1];
rz(-0.70530546) q[1];
sx q[1];
rz(-0.7440872) q[1];
rz(-0.75833851) q[3];
sx q[3];
rz(-1.3876649) q[3];
sx q[3];
rz(1.5116215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4522176) q[2];
sx q[2];
rz(-1.3000501) q[2];
sx q[2];
rz(-2.8049862) q[2];
rz(1.5161139) q[3];
sx q[3];
rz(-2.5879526) q[3];
sx q[3];
rz(1.5256933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7933554) q[0];
sx q[0];
rz(-2.0331148) q[0];
sx q[0];
rz(-0.021214699) q[0];
rz(1.9477828) q[1];
sx q[1];
rz(-2.1021011) q[1];
sx q[1];
rz(-2.3056727) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7556691) q[0];
sx q[0];
rz(-0.14005157) q[0];
sx q[0];
rz(1.7007909) q[0];
rz(-pi) q[1];
rz(2.1875728) q[2];
sx q[2];
rz(-2.4980133) q[2];
sx q[2];
rz(1.8804903) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1509717) q[1];
sx q[1];
rz(-2.4358106) q[1];
sx q[1];
rz(1.1433931) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7153035) q[3];
sx q[3];
rz(-2.3549035) q[3];
sx q[3];
rz(-2.8601437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8643643) q[2];
sx q[2];
rz(-1.9936864) q[2];
sx q[2];
rz(-1.345984) q[2];
rz(-2.7820382) q[3];
sx q[3];
rz(-0.94272009) q[3];
sx q[3];
rz(-2.6446222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42831746) q[0];
sx q[0];
rz(-2.064216) q[0];
sx q[0];
rz(-2.0879478) q[0];
rz(1.2288278) q[1];
sx q[1];
rz(-1.6002974) q[1];
sx q[1];
rz(0.4371117) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29262221) q[0];
sx q[0];
rz(-1.4842352) q[0];
sx q[0];
rz(-1.863088) q[0];
rz(-pi) q[1];
rz(-2.897981) q[2];
sx q[2];
rz(-1.1738452) q[2];
sx q[2];
rz(-1.5872019) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.25698173) q[1];
sx q[1];
rz(-1.4323438) q[1];
sx q[1];
rz(2.715766) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4594853) q[3];
sx q[3];
rz(-2.6442332) q[3];
sx q[3];
rz(1.5745844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.019471021) q[2];
sx q[2];
rz(-0.78142587) q[2];
sx q[2];
rz(-1.0220698) q[2];
rz(1.2381037) q[3];
sx q[3];
rz(-0.3823897) q[3];
sx q[3];
rz(-0.4195655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5220752) q[0];
sx q[0];
rz(-1.2457122) q[0];
sx q[0];
rz(-2.1602901) q[0];
rz(-3.006382) q[1];
sx q[1];
rz(-2.0573261) q[1];
sx q[1];
rz(0.19128004) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6387588) q[0];
sx q[0];
rz(-1.6573818) q[0];
sx q[0];
rz(1.4214574) q[0];
rz(-pi) q[1];
rz(-2.3046266) q[2];
sx q[2];
rz(-1.3651197) q[2];
sx q[2];
rz(-0.93081805) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7820218) q[1];
sx q[1];
rz(-1.4683717) q[1];
sx q[1];
rz(-0.77918474) q[1];
x q[2];
rz(-1.220827) q[3];
sx q[3];
rz(-1.063949) q[3];
sx q[3];
rz(0.99196539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.68025756) q[2];
sx q[2];
rz(-0.98538435) q[2];
sx q[2];
rz(-2.130924) q[2];
rz(-2.3800395) q[3];
sx q[3];
rz(-1.1798309) q[3];
sx q[3];
rz(0.23553577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8028832) q[0];
sx q[0];
rz(-2.8864679) q[0];
sx q[0];
rz(-0.55661911) q[0];
rz(3.026475) q[1];
sx q[1];
rz(-1.8042253) q[1];
sx q[1];
rz(-2.1690878) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2199729) q[0];
sx q[0];
rz(-3.0191506) q[0];
sx q[0];
rz(0.58453154) q[0];
rz(0.33072492) q[2];
sx q[2];
rz(-2.2933368) q[2];
sx q[2];
rz(-1.5028138) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.433686) q[1];
sx q[1];
rz(-1.7644595) q[1];
sx q[1];
rz(0.14240264) q[1];
rz(0.036690849) q[3];
sx q[3];
rz(-2.1864236) q[3];
sx q[3];
rz(0.43945593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.82289034) q[2];
sx q[2];
rz(-1.0914785) q[2];
sx q[2];
rz(1.5931607) q[2];
rz(1.3657773) q[3];
sx q[3];
rz(-0.32309353) q[3];
sx q[3];
rz(0.95388609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71972972) q[0];
sx q[0];
rz(-1.8780163) q[0];
sx q[0];
rz(1.7156037) q[0];
rz(1.0643719) q[1];
sx q[1];
rz(-1.0168889) q[1];
sx q[1];
rz(2.7672966) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9126041) q[0];
sx q[0];
rz(-1.2545663) q[0];
sx q[0];
rz(2.7094748) q[0];
x q[1];
rz(2.6475545) q[2];
sx q[2];
rz(-1.8015773) q[2];
sx q[2];
rz(-1.5649232) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4243851) q[1];
sx q[1];
rz(-1.8339515) q[1];
sx q[1];
rz(-0.56916635) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3338942) q[3];
sx q[3];
rz(-2.2904615) q[3];
sx q[3];
rz(2.6157275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2465683) q[2];
sx q[2];
rz(-2.4763069) q[2];
sx q[2];
rz(0.95823112) q[2];
rz(0.22917497) q[3];
sx q[3];
rz(-1.4567679) q[3];
sx q[3];
rz(-2.5206101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7763057) q[0];
sx q[0];
rz(-1.1927274) q[0];
sx q[0];
rz(2.2348485) q[0];
rz(-1.0892185) q[1];
sx q[1];
rz(-1.4995432) q[1];
sx q[1];
rz(1.3100756) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.037127) q[0];
sx q[0];
rz(-1.6220399) q[0];
sx q[0];
rz(-0.38802223) q[0];
rz(-0.71307619) q[2];
sx q[2];
rz(-1.3165858) q[2];
sx q[2];
rz(-1.6377246) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.20565614) q[1];
sx q[1];
rz(-1.2381136) q[1];
sx q[1];
rz(1.7722305) q[1];
rz(1.5982315) q[3];
sx q[3];
rz(-2.1260288) q[3];
sx q[3];
rz(-2.1813986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92581302) q[2];
sx q[2];
rz(-0.44162193) q[2];
sx q[2];
rz(-1.4833935) q[2];
rz(-0.27967927) q[3];
sx q[3];
rz(-2.1488583) q[3];
sx q[3];
rz(0.057597615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9782372) q[0];
sx q[0];
rz(-1.5196479) q[0];
sx q[0];
rz(-0.21959198) q[0];
rz(-0.50312463) q[1];
sx q[1];
rz(-2.2527835) q[1];
sx q[1];
rz(-0.84987744) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4655315) q[0];
sx q[0];
rz(-1.413835) q[0];
sx q[0];
rz(-1.3423052) q[0];
rz(-pi) q[1];
rz(1.6860028) q[2];
sx q[2];
rz(-0.68636471) q[2];
sx q[2];
rz(-2.3442868) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1666607) q[1];
sx q[1];
rz(-0.76172511) q[1];
sx q[1];
rz(-1.5207661) q[1];
rz(-pi) q[2];
rz(1.1342808) q[3];
sx q[3];
rz(-2.8562162) q[3];
sx q[3];
rz(-0.34261045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5510817) q[2];
sx q[2];
rz(-1.8293646) q[2];
sx q[2];
rz(-1.3809416) q[2];
rz(0.75602174) q[3];
sx q[3];
rz(-2.9383926) q[3];
sx q[3];
rz(2.7856564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6324156) q[0];
sx q[0];
rz(-2.248705) q[0];
sx q[0];
rz(0.40503043) q[0];
rz(-2.6889154) q[1];
sx q[1];
rz(-2.15937) q[1];
sx q[1];
rz(1.2776432) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1179827) q[0];
sx q[0];
rz(-1.3408459) q[0];
sx q[0];
rz(0.50595565) q[0];
rz(-2.5759376) q[2];
sx q[2];
rz(-2.3662162) q[2];
sx q[2];
rz(-2.2765809) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.35754044) q[1];
sx q[1];
rz(-1.6077542) q[1];
sx q[1];
rz(1.1251015) q[1];
rz(-pi) q[2];
rz(-1.9926662) q[3];
sx q[3];
rz(-2.359169) q[3];
sx q[3];
rz(-0.22280927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8032802) q[2];
sx q[2];
rz(-2.7351604) q[2];
sx q[2];
rz(2.7837616) q[2];
rz(1.4194277) q[3];
sx q[3];
rz(-1.2737041) q[3];
sx q[3];
rz(-2.0675802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2232067) q[0];
sx q[0];
rz(-0.077843852) q[0];
sx q[0];
rz(0.11225587) q[0];
rz(-0.90011251) q[1];
sx q[1];
rz(-1.0670412) q[1];
sx q[1];
rz(0.21044883) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99104098) q[0];
sx q[0];
rz(-1.3897087) q[0];
sx q[0];
rz(0.51245706) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9774882) q[2];
sx q[2];
rz(-1.9774984) q[2];
sx q[2];
rz(-1.1526398) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.12293832) q[1];
sx q[1];
rz(-1.2631053) q[1];
sx q[1];
rz(1.7880746) q[1];
x q[2];
rz(-1.3335161) q[3];
sx q[3];
rz(-1.8144326) q[3];
sx q[3];
rz(-0.39792774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9615053) q[2];
sx q[2];
rz(-2.4752361) q[2];
sx q[2];
rz(1.5562742) q[2];
rz(-1.2735584) q[3];
sx q[3];
rz(-2.5189416) q[3];
sx q[3];
rz(-2.9343228) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56959854) q[0];
sx q[0];
rz(-0.8710237) q[0];
sx q[0];
rz(-1.3652753) q[0];
rz(-0.81644425) q[1];
sx q[1];
rz(-1.2533617) q[1];
sx q[1];
rz(-0.15773699) q[1];
rz(1.1007166) q[2];
sx q[2];
rz(-1.7190949) q[2];
sx q[2];
rz(0.20863056) q[2];
rz(1.9498701) q[3];
sx q[3];
rz(-2.3754397) q[3];
sx q[3];
rz(0.55879186) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
