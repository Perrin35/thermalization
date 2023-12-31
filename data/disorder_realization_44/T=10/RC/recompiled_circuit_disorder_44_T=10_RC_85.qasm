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
rz(6.7232806) q[0];
sx q[0];
rz(6.4203782) q[0];
rz(1.4057012) q[1];
sx q[1];
rz(4.5448137) q[1];
sx q[1];
rz(9.9546976) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0041381) q[0];
sx q[0];
rz(-1.9507427) q[0];
sx q[0];
rz(3.0259892) q[0];
rz(-pi) q[1];
rz(-0.86962236) q[2];
sx q[2];
rz(-2.6343971) q[2];
sx q[2];
rz(1.6091572) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.73218988) q[1];
sx q[1];
rz(-0.70530546) q[1];
sx q[1];
rz(2.3975055) q[1];
x q[2];
rz(-2.8785273) q[3];
sx q[3];
rz(-2.3657626) q[3];
sx q[3];
rz(2.8924931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.68937504) q[2];
sx q[2];
rz(-1.3000501) q[2];
sx q[2];
rz(-2.8049862) q[2];
rz(-1.6254788) q[3];
sx q[3];
rz(-2.5879526) q[3];
sx q[3];
rz(-1.6158993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7933554) q[0];
sx q[0];
rz(-2.0331148) q[0];
sx q[0];
rz(3.120378) q[0];
rz(-1.9477828) q[1];
sx q[1];
rz(-1.0394916) q[1];
sx q[1];
rz(-2.3056727) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5171889) q[0];
sx q[0];
rz(-1.7096585) q[0];
sx q[0];
rz(-3.1233203) q[0];
x q[1];
rz(-2.1199273) q[2];
sx q[2];
rz(-1.9252535) q[2];
sx q[2];
rz(-0.82565386) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6902496) q[1];
sx q[1];
rz(-2.2022044) q[1];
sx q[1];
rz(0.33957014) q[1];
rz(-pi) q[2];
rz(-0.42628916) q[3];
sx q[3];
rz(-2.3549035) q[3];
sx q[3];
rz(-2.8601437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8643643) q[2];
sx q[2];
rz(-1.1479062) q[2];
sx q[2];
rz(-1.7956087) q[2];
rz(-2.7820382) q[3];
sx q[3];
rz(-2.1988726) q[3];
sx q[3];
rz(-0.49697044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.7132752) q[0];
sx q[0];
rz(-2.064216) q[0];
sx q[0];
rz(-1.0536449) q[0];
rz(1.9127649) q[1];
sx q[1];
rz(-1.6002974) q[1];
sx q[1];
rz(-0.4371117) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29262221) q[0];
sx q[0];
rz(-1.4842352) q[0];
sx q[0];
rz(-1.863088) q[0];
rz(-pi) q[1];
x q[1];
rz(1.978546) q[2];
sx q[2];
rz(-1.3464658) q[2];
sx q[2];
rz(-3.0293904) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6091842) q[1];
sx q[1];
rz(-0.44645616) q[1];
sx q[1];
rz(-2.8162454) q[1];
x q[2];
rz(1.4594853) q[3];
sx q[3];
rz(-0.49735945) q[3];
sx q[3];
rz(-1.5745844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1221216) q[2];
sx q[2];
rz(-0.78142587) q[2];
sx q[2];
rz(2.1195228) q[2];
rz(1.9034889) q[3];
sx q[3];
rz(-2.759203) q[3];
sx q[3];
rz(-0.4195655) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5220752) q[0];
sx q[0];
rz(-1.2457122) q[0];
sx q[0];
rz(0.98130256) q[0];
rz(-0.13521067) q[1];
sx q[1];
rz(-1.0842666) q[1];
sx q[1];
rz(-2.9503126) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.080973074) q[0];
sx q[0];
rz(-1.7195716) q[0];
sx q[0];
rz(0.087555126) q[0];
rz(-2.8677167) q[2];
sx q[2];
rz(-0.855815) q[2];
sx q[2];
rz(-2.3194734) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0335238) q[1];
sx q[1];
rz(-2.3571157) q[1];
sx q[1];
rz(0.14524059) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6077765) q[3];
sx q[3];
rz(-1.2663519) q[3];
sx q[3];
rz(0.75418762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.68025756) q[2];
sx q[2];
rz(-2.1562083) q[2];
sx q[2];
rz(1.0106687) q[2];
rz(2.3800395) q[3];
sx q[3];
rz(-1.1798309) q[3];
sx q[3];
rz(2.9060569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33870944) q[0];
sx q[0];
rz(-2.8864679) q[0];
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
rz(0.92161979) q[0];
sx q[0];
rz(-0.12244206) q[0];
sx q[0];
rz(-2.5570611) q[0];
x q[1];
rz(2.8108677) q[2];
sx q[2];
rz(-0.84825584) q[2];
sx q[2];
rz(1.6387788) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7079066) q[1];
sx q[1];
rz(-1.7644595) q[1];
sx q[1];
rz(-2.99919) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6226107) q[3];
sx q[3];
rz(-0.6165781) q[3];
sx q[3];
rz(-0.50293621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82289034) q[2];
sx q[2];
rz(-1.0914785) q[2];
sx q[2];
rz(-1.548432) q[2];
rz(1.7758153) q[3];
sx q[3];
rz(-2.8184991) q[3];
sx q[3];
rz(0.95388609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4218629) q[0];
sx q[0];
rz(-1.8780163) q[0];
sx q[0];
rz(1.7156037) q[0];
rz(2.0772207) q[1];
sx q[1];
rz(-1.0168889) q[1];
sx q[1];
rz(-2.7672966) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4842589) q[0];
sx q[0];
rz(-1.1614292) q[0];
sx q[0];
rz(-1.9166458) q[0];
x q[1];
rz(0.49403814) q[2];
sx q[2];
rz(-1.3400153) q[2];
sx q[2];
rz(1.5766694) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.018505521) q[1];
sx q[1];
rz(-1.0235041) q[1];
sx q[1];
rz(1.8803384) q[1];
x q[2];
rz(2.3338942) q[3];
sx q[3];
rz(-0.85113111) q[3];
sx q[3];
rz(0.52586517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8950243) q[2];
sx q[2];
rz(-2.4763069) q[2];
sx q[2];
rz(2.1833615) q[2];
rz(-0.22917497) q[3];
sx q[3];
rz(-1.4567679) q[3];
sx q[3];
rz(-0.62098256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7763057) q[0];
sx q[0];
rz(-1.1927274) q[0];
sx q[0];
rz(-2.2348485) q[0];
rz(1.0892185) q[1];
sx q[1];
rz(-1.4995432) q[1];
sx q[1];
rz(-1.3100756) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7999254) q[0];
sx q[0];
rz(-2.7503715) q[0];
sx q[0];
rz(3.0068586) q[0];
rz(-0.37808772) q[2];
sx q[2];
rz(-2.3921161) q[2];
sx q[2];
rz(-2.7917002) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3771364) q[1];
sx q[1];
rz(-0.38696445) q[1];
sx q[1];
rz(-0.52486921) q[1];
rz(-pi) q[2];
rz(-1.5433611) q[3];
sx q[3];
rz(-1.0155639) q[3];
sx q[3];
rz(-0.96019402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.92581302) q[2];
sx q[2];
rz(-2.6999707) q[2];
sx q[2];
rz(1.4833935) q[2];
rz(0.27967927) q[3];
sx q[3];
rz(-0.99273434) q[3];
sx q[3];
rz(-3.083995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16335547) q[0];
sx q[0];
rz(-1.5196479) q[0];
sx q[0];
rz(-0.21959198) q[0];
rz(-2.638468) q[1];
sx q[1];
rz(-2.2527835) q[1];
sx q[1];
rz(0.84987744) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6760611) q[0];
sx q[0];
rz(-1.413835) q[0];
sx q[0];
rz(-1.7992875) q[0];
x q[1];
rz(-3.0476961) q[2];
sx q[2];
rz(-2.2517423) q[2];
sx q[2];
rz(0.64881334) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.5596458) q[1];
sx q[1];
rz(-1.5362745) q[1];
sx q[1];
rz(2.3318961) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12340801) q[3];
sx q[3];
rz(-1.828769) q[3];
sx q[3];
rz(-0.79515275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5510817) q[2];
sx q[2];
rz(-1.8293646) q[2];
sx q[2];
rz(1.3809416) q[2];
rz(0.75602174) q[3];
sx q[3];
rz(-2.9383926) q[3];
sx q[3];
rz(-0.35593629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6324156) q[0];
sx q[0];
rz(-2.248705) q[0];
sx q[0];
rz(2.7365622) q[0];
rz(-2.6889154) q[1];
sx q[1];
rz(-2.15937) q[1];
sx q[1];
rz(1.2776432) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1179827) q[0];
sx q[0];
rz(-1.8007468) q[0];
sx q[0];
rz(0.50595565) q[0];
rz(-pi) q[1];
rz(-2.4503166) q[2];
sx q[2];
rz(-1.1862159) q[2];
sx q[2];
rz(-2.861475) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7840522) q[1];
sx q[1];
rz(-1.5338384) q[1];
sx q[1];
rz(1.1251015) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9926662) q[3];
sx q[3];
rz(-0.78242362) q[3];
sx q[3];
rz(2.9187834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8032802) q[2];
sx q[2];
rz(-2.7351604) q[2];
sx q[2];
rz(-0.35783106) q[2];
rz(-1.7221649) q[3];
sx q[3];
rz(-1.8678886) q[3];
sx q[3];
rz(-1.0740124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.91838592) q[0];
sx q[0];
rz(-3.0637488) q[0];
sx q[0];
rz(0.11225587) q[0];
rz(0.90011251) q[1];
sx q[1];
rz(-1.0670412) q[1];
sx q[1];
rz(-0.21044883) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2699379) q[0];
sx q[0];
rz(-2.6007814) q[0];
sx q[0];
rz(0.35738118) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9774882) q[2];
sx q[2];
rz(-1.1640942) q[2];
sx q[2];
rz(1.1526398) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0186543) q[1];
sx q[1];
rz(-1.8784874) q[1];
sx q[1];
rz(1.3535181) q[1];
x q[2];
rz(-2.8912192) q[3];
sx q[3];
rz(-1.3406521) q[3];
sx q[3];
rz(1.2311414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9615053) q[2];
sx q[2];
rz(-0.6663565) q[2];
sx q[2];
rz(-1.5853184) q[2];
rz(1.8680343) q[3];
sx q[3];
rz(-2.5189416) q[3];
sx q[3];
rz(-2.9343228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5719941) q[0];
sx q[0];
rz(-2.270569) q[0];
sx q[0];
rz(1.7763174) q[0];
rz(2.3251484) q[1];
sx q[1];
rz(-1.2533617) q[1];
sx q[1];
rz(-0.15773699) q[1];
rz(1.1007166) q[2];
sx q[2];
rz(-1.7190949) q[2];
sx q[2];
rz(0.20863056) q[2];
rz(0.84135009) q[3];
sx q[3];
rz(-1.311306) q[3];
sx q[3];
rz(-0.73248274) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
