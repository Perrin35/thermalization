OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8496348) q[0];
sx q[0];
rz(-0.44009527) q[0];
sx q[0];
rz(0.13719288) q[0];
rz(1.4057012) q[1];
sx q[1];
rz(-1.7383716) q[1];
sx q[1];
rz(0.52991968) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44088988) q[0];
sx q[0];
rz(-0.39632495) q[0];
sx q[0];
rz(-1.8519782) q[0];
rz(-pi) q[1];
rz(2.797384) q[2];
sx q[2];
rz(-1.95103) q[2];
sx q[2];
rz(-0.76438475) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5296386) q[1];
sx q[1];
rz(-1.0736335) q[1];
sx q[1];
rz(-1.0477209) q[1];
rz(-2.3832541) q[3];
sx q[3];
rz(-1.7539277) q[3];
sx q[3];
rz(-1.6299712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4522176) q[2];
sx q[2];
rz(-1.3000501) q[2];
sx q[2];
rz(2.8049862) q[2];
rz(1.5161139) q[3];
sx q[3];
rz(-0.55364004) q[3];
sx q[3];
rz(-1.5256933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7933554) q[0];
sx q[0];
rz(-1.1084778) q[0];
sx q[0];
rz(0.021214699) q[0];
rz(-1.9477828) q[1];
sx q[1];
rz(-1.0394916) q[1];
sx q[1];
rz(-2.3056727) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5171889) q[0];
sx q[0];
rz(-1.7096585) q[0];
sx q[0];
rz(-0.018272321) q[0];
rz(-pi) q[1];
rz(-1.0216653) q[2];
sx q[2];
rz(-1.2163391) q[2];
sx q[2];
rz(-0.82565386) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1509717) q[1];
sx q[1];
rz(-2.4358106) q[1];
sx q[1];
rz(-1.1433931) q[1];
x q[2];
rz(2.4017176) q[3];
sx q[3];
rz(-1.8679108) q[3];
sx q[3];
rz(1.5418996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2772284) q[2];
sx q[2];
rz(-1.1479062) q[2];
sx q[2];
rz(-1.345984) q[2];
rz(-2.7820382) q[3];
sx q[3];
rz(-2.1988726) q[3];
sx q[3];
rz(-0.49697044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7132752) q[0];
sx q[0];
rz(-2.064216) q[0];
sx q[0];
rz(-2.0879478) q[0];
rz(1.2288278) q[1];
sx q[1];
rz(-1.6002974) q[1];
sx q[1];
rz(0.4371117) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8374098) q[0];
sx q[0];
rz(-1.8619616) q[0];
sx q[0];
rz(3.0512179) q[0];
rz(-pi) q[1];
rz(2.0929298) q[2];
sx q[2];
rz(-0.46233593) q[2];
sx q[2];
rz(-0.98302746) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8902953) q[1];
sx q[1];
rz(-1.9922868) q[1];
sx q[1];
rz(1.7226268) q[1];
rz(-pi) q[2];
rz(-2.0655572) q[3];
sx q[3];
rz(-1.5177739) q[3];
sx q[3];
rz(-3.0474636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1221216) q[2];
sx q[2];
rz(-2.3601668) q[2];
sx q[2];
rz(-1.0220698) q[2];
rz(1.9034889) q[3];
sx q[3];
rz(-0.3823897) q[3];
sx q[3];
rz(-2.7220272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5220752) q[0];
sx q[0];
rz(-1.2457122) q[0];
sx q[0];
rz(2.1602901) q[0];
rz(3.006382) q[1];
sx q[1];
rz(-1.0842666) q[1];
sx q[1];
rz(-2.9503126) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0606196) q[0];
sx q[0];
rz(-1.7195716) q[0];
sx q[0];
rz(3.0540375) q[0];
x q[1];
rz(0.27387597) q[2];
sx q[2];
rz(-0.855815) q[2];
sx q[2];
rz(-2.3194734) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3595708) q[1];
sx q[1];
rz(-1.4683717) q[1];
sx q[1];
rz(0.77918474) q[1];
x q[2];
rz(2.6077765) q[3];
sx q[3];
rz(-1.2663519) q[3];
sx q[3];
rz(2.387405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4613351) q[2];
sx q[2];
rz(-0.98538435) q[2];
sx q[2];
rz(1.0106687) q[2];
rz(-2.3800395) q[3];
sx q[3];
rz(-1.1798309) q[3];
sx q[3];
rz(0.23553577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8028832) q[0];
sx q[0];
rz(-0.25512472) q[0];
sx q[0];
rz(2.5849735) q[0];
rz(0.11511766) q[1];
sx q[1];
rz(-1.8042253) q[1];
sx q[1];
rz(-0.97250485) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33362493) q[0];
sx q[0];
rz(-1.4687612) q[0];
sx q[0];
rz(1.5029961) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9237256) q[2];
sx q[2];
rz(-2.3595516) q[2];
sx q[2];
rz(-2.1176586) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.16469615) q[1];
sx q[1];
rz(-1.7105192) q[1];
sx q[1];
rz(1.7663899) q[1];
x q[2];
rz(-1.5189819) q[3];
sx q[3];
rz(-0.6165781) q[3];
sx q[3];
rz(-2.6386564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3187023) q[2];
sx q[2];
rz(-2.0501142) q[2];
sx q[2];
rz(1.548432) q[2];
rz(-1.7758153) q[3];
sx q[3];
rz(-0.32309353) q[3];
sx q[3];
rz(0.95388609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(2.4218629) q[0];
sx q[0];
rz(-1.8780163) q[0];
sx q[0];
rz(-1.7156037) q[0];
rz(1.0643719) q[1];
sx q[1];
rz(-2.1247037) q[1];
sx q[1];
rz(-2.7672966) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4842589) q[0];
sx q[0];
rz(-1.1614292) q[0];
sx q[0];
rz(1.9166458) q[0];
rz(-2.6815368) q[2];
sx q[2];
rz(-2.6003777) q[2];
sx q[2];
rz(-0.40749007) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4243851) q[1];
sx q[1];
rz(-1.8339515) q[1];
sx q[1];
rz(-2.5724263) q[1];
rz(-pi) q[2];
rz(-2.3338942) q[3];
sx q[3];
rz(-0.85113111) q[3];
sx q[3];
rz(2.6157275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2465683) q[2];
sx q[2];
rz(-2.4763069) q[2];
sx q[2];
rz(-2.1833615) q[2];
rz(0.22917497) q[3];
sx q[3];
rz(-1.6848247) q[3];
sx q[3];
rz(-0.62098256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7763057) q[0];
sx q[0];
rz(-1.9488652) q[0];
sx q[0];
rz(-2.2348485) q[0];
rz(1.0892185) q[1];
sx q[1];
rz(-1.4995432) q[1];
sx q[1];
rz(-1.3100756) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7999254) q[0];
sx q[0];
rz(-0.39122117) q[0];
sx q[0];
rz(0.13473405) q[0];
rz(0.37808772) q[2];
sx q[2];
rz(-2.3921161) q[2];
sx q[2];
rz(2.7917002) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.20565614) q[1];
sx q[1];
rz(-1.2381136) q[1];
sx q[1];
rz(1.3693621) q[1];
rz(-pi) q[2];
rz(-0.044192627) q[3];
sx q[3];
rz(-2.5857539) q[3];
sx q[3];
rz(1.0122055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2157796) q[2];
sx q[2];
rz(-0.44162193) q[2];
sx q[2];
rz(-1.4833935) q[2];
rz(-0.27967927) q[3];
sx q[3];
rz(-0.99273434) q[3];
sx q[3];
rz(3.083995) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
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
rz(-0.88880912) q[1];
sx q[1];
rz(-2.2917152) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4655315) q[0];
sx q[0];
rz(-1.7277576) q[0];
sx q[0];
rz(1.7992875) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4555898) q[2];
sx q[2];
rz(-2.4552279) q[2];
sx q[2];
rz(-0.79730588) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0440158) q[1];
sx q[1];
rz(-0.8102639) q[1];
sx q[1];
rz(-3.0939328) q[1];
rz(-pi) q[2];
rz(-0.12340801) q[3];
sx q[3];
rz(-1.828769) q[3];
sx q[3];
rz(-2.3464399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5510817) q[2];
sx q[2];
rz(-1.8293646) q[2];
sx q[2];
rz(1.760651) q[2];
rz(-0.75602174) q[3];
sx q[3];
rz(-2.9383926) q[3];
sx q[3];
rz(-2.7856564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5091771) q[0];
sx q[0];
rz(-2.248705) q[0];
sx q[0];
rz(-2.7365622) q[0];
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
rz(1.02361) q[0];
sx q[0];
rz(-1.3408459) q[0];
sx q[0];
rz(0.50595565) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56565506) q[2];
sx q[2];
rz(-2.3662162) q[2];
sx q[2];
rz(2.2765809) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2904418) q[1];
sx q[1];
rz(-2.6944707) q[1];
sx q[1];
rz(1.4852344) q[1];
x q[2];
rz(1.9926662) q[3];
sx q[3];
rz(-0.78242362) q[3];
sx q[3];
rz(-0.22280927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3383125) q[2];
sx q[2];
rz(-0.40643224) q[2];
sx q[2];
rz(-2.7837616) q[2];
rz(1.4194277) q[3];
sx q[3];
rz(-1.8678886) q[3];
sx q[3];
rz(2.0675802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91838592) q[0];
sx q[0];
rz(-3.0637488) q[0];
sx q[0];
rz(3.0293368) q[0];
rz(-2.2414801) q[1];
sx q[1];
rz(-2.0745514) q[1];
sx q[1];
rz(0.21044883) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68073273) q[0];
sx q[0];
rz(-1.0675149) q[0];
sx q[0];
rz(-1.3637278) q[0];
rz(2.9774882) q[2];
sx q[2];
rz(-1.9774984) q[2];
sx q[2];
rz(1.1526398) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.12293832) q[1];
sx q[1];
rz(-1.8784874) q[1];
sx q[1];
rz(-1.7880746) q[1];
rz(-pi) q[2];
rz(-2.8912192) q[3];
sx q[3];
rz(-1.3406521) q[3];
sx q[3];
rz(1.2311414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9615053) q[2];
sx q[2];
rz(-0.6663565) q[2];
sx q[2];
rz(-1.5562742) q[2];
rz(1.8680343) q[3];
sx q[3];
rz(-0.62265101) q[3];
sx q[3];
rz(-0.20726985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56959854) q[0];
sx q[0];
rz(-2.270569) q[0];
sx q[0];
rz(1.7763174) q[0];
rz(-0.81644425) q[1];
sx q[1];
rz(-1.2533617) q[1];
sx q[1];
rz(-0.15773699) q[1];
rz(2.9755637) q[2];
sx q[2];
rz(-1.1062853) q[2];
sx q[2];
rz(1.7044978) q[2];
rz(2.3002426) q[3];
sx q[3];
rz(-1.8302866) q[3];
sx q[3];
rz(2.4091099) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];