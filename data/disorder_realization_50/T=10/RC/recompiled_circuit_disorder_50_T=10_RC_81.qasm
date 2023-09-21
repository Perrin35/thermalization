OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.87941909) q[0];
sx q[0];
rz(4.8882422) q[0];
sx q[0];
rz(12.56765) q[0];
rz(-1.6969504) q[1];
sx q[1];
rz(-2.0386219) q[1];
sx q[1];
rz(-2.3666518) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9026731) q[0];
sx q[0];
rz(-2.9905149) q[0];
sx q[0];
rz(-2.8047049) q[0];
rz(-pi) q[1];
rz(-2.7156419) q[2];
sx q[2];
rz(-0.89086878) q[2];
sx q[2];
rz(-1.8217063) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.70667627) q[1];
sx q[1];
rz(-1.4286563) q[1];
sx q[1];
rz(3.1013156) q[1];
x q[2];
rz(-1.7232056) q[3];
sx q[3];
rz(-1.8525575) q[3];
sx q[3];
rz(-1.2581274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0455735) q[2];
sx q[2];
rz(-1.1117659) q[2];
sx q[2];
rz(-1.1958896) q[2];
rz(-1.1536417) q[3];
sx q[3];
rz(-0.78912815) q[3];
sx q[3];
rz(1.7529863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7213223) q[0];
sx q[0];
rz(-0.35636154) q[0];
sx q[0];
rz(-1.2715682) q[0];
rz(-1.0999854) q[1];
sx q[1];
rz(-1.1106691) q[1];
sx q[1];
rz(-1.7659448) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9346416) q[0];
sx q[0];
rz(-1.9793545) q[0];
sx q[0];
rz(-1.4757399) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8091082) q[2];
sx q[2];
rz(-1.6147699) q[2];
sx q[2];
rz(-2.6251052) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.34054204) q[1];
sx q[1];
rz(-1.2060296) q[1];
sx q[1];
rz(1.3786475) q[1];
rz(1.0166753) q[3];
sx q[3];
rz(-1.2614054) q[3];
sx q[3];
rz(2.4966937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9282844) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(-2.6518872) q[2];
rz(-2.0080163) q[3];
sx q[3];
rz(-2.9462892) q[3];
sx q[3];
rz(1.0666696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5415444) q[0];
sx q[0];
rz(-1.2788037) q[0];
sx q[0];
rz(-2.9009853) q[0];
rz(2.799017) q[1];
sx q[1];
rz(-2.1668285) q[1];
sx q[1];
rz(-1.2352357) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2152527) q[0];
sx q[0];
rz(-1.0010166) q[0];
sx q[0];
rz(2.0195872) q[0];
rz(-pi) q[1];
rz(-2.5580102) q[2];
sx q[2];
rz(-0.43791134) q[2];
sx q[2];
rz(-2.8180518) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0864799) q[1];
sx q[1];
rz(-1.8243454) q[1];
sx q[1];
rz(-2.8716645) q[1];
rz(-pi) q[2];
rz(-2.5123185) q[3];
sx q[3];
rz(-2.4089703) q[3];
sx q[3];
rz(-0.34793684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.76413313) q[2];
sx q[2];
rz(-1.5185792) q[2];
sx q[2];
rz(1.7049449) q[2];
rz(1.7403587) q[3];
sx q[3];
rz(-1.2652206) q[3];
sx q[3];
rz(0.60825545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075994611) q[0];
sx q[0];
rz(-1.6125212) q[0];
sx q[0];
rz(-2.3377989) q[0];
rz(-2.1919788) q[1];
sx q[1];
rz(-1.6751553) q[1];
sx q[1];
rz(-3.0217357) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6896497) q[0];
sx q[0];
rz(-1.6065238) q[0];
sx q[0];
rz(-2.1769051) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7718133) q[2];
sx q[2];
rz(-0.38447194) q[2];
sx q[2];
rz(-2.0535) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3164697) q[1];
sx q[1];
rz(-0.23288647) q[1];
sx q[1];
rz(-1.6474001) q[1];
x q[2];
rz(2.5424764) q[3];
sx q[3];
rz(-1.1154419) q[3];
sx q[3];
rz(-1.276254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5084761) q[2];
sx q[2];
rz(-1.2458331) q[2];
sx q[2];
rz(-0.072337739) q[2];
rz(-2.7667601) q[3];
sx q[3];
rz(-0.65632498) q[3];
sx q[3];
rz(1.9434631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4478093) q[0];
sx q[0];
rz(-0.57149514) q[0];
sx q[0];
rz(2.6547292) q[0];
rz(2.411719) q[1];
sx q[1];
rz(-0.90881538) q[1];
sx q[1];
rz(1.240085) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86868984) q[0];
sx q[0];
rz(-1.5839974) q[0];
sx q[0];
rz(-0.23916434) q[0];
x q[1];
rz(-1.8334332) q[2];
sx q[2];
rz(-1.8034168) q[2];
sx q[2];
rz(0.55839415) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1357437) q[1];
sx q[1];
rz(-1.504335) q[1];
sx q[1];
rz(2.5367516) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2206743) q[3];
sx q[3];
rz(-0.71483597) q[3];
sx q[3];
rz(-2.9361847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.038625) q[2];
sx q[2];
rz(-0.90927783) q[2];
sx q[2];
rz(-0.70303482) q[2];
rz(1.7317584) q[3];
sx q[3];
rz(-1.3491646) q[3];
sx q[3];
rz(-2.9838802) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96419656) q[0];
sx q[0];
rz(-1.2921768) q[0];
sx q[0];
rz(-2.1512206) q[0];
rz(3.0888427) q[1];
sx q[1];
rz(-2.2149448) q[1];
sx q[1];
rz(1.4809158) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4617417) q[0];
sx q[0];
rz(-1.4977507) q[0];
sx q[0];
rz(-1.1887656) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5691721) q[2];
sx q[2];
rz(-0.44518984) q[2];
sx q[2];
rz(-2.0507398) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.33038352) q[1];
sx q[1];
rz(-0.57302176) q[1];
sx q[1];
rz(-1.9622383) q[1];
rz(3.092993) q[3];
sx q[3];
rz(-2.2489293) q[3];
sx q[3];
rz(2.5130659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.0084373077) q[2];
sx q[2];
rz(-1.4921654) q[2];
sx q[2];
rz(-0.56224242) q[2];
rz(-1.0605313) q[3];
sx q[3];
rz(-0.74354592) q[3];
sx q[3];
rz(-0.26091584) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9437207) q[0];
sx q[0];
rz(-1.7892388) q[0];
sx q[0];
rz(1.6725756) q[0];
rz(-2.127227) q[1];
sx q[1];
rz(-1.0275774) q[1];
sx q[1];
rz(1.8168824) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4949188) q[0];
sx q[0];
rz(-3.0122628) q[0];
sx q[0];
rz(-2.9652251) q[0];
x q[1];
rz(-1.0075188) q[2];
sx q[2];
rz(-2.0118606) q[2];
sx q[2];
rz(1.2457459) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.087254698) q[1];
sx q[1];
rz(-1.9458276) q[1];
sx q[1];
rz(0.43941984) q[1];
rz(-1.0876098) q[3];
sx q[3];
rz(-2.4251068) q[3];
sx q[3];
rz(-2.0487752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5801195) q[2];
sx q[2];
rz(-1.3886398) q[2];
sx q[2];
rz(-0.44357792) q[2];
rz(0.94868547) q[3];
sx q[3];
rz(-1.1921927) q[3];
sx q[3];
rz(2.366812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7396486) q[0];
sx q[0];
rz(-1.0111324) q[0];
sx q[0];
rz(0.46052128) q[0];
rz(3.0415688) q[1];
sx q[1];
rz(-2.1477551) q[1];
sx q[1];
rz(-1.8519648) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78532082) q[0];
sx q[0];
rz(-0.7542146) q[0];
sx q[0];
rz(0.69099364) q[0];
rz(-2.4256267) q[2];
sx q[2];
rz(-2.1858474) q[2];
sx q[2];
rz(-1.3032608) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.21266567) q[1];
sx q[1];
rz(-1.3839098) q[1];
sx q[1];
rz(-3.0287663) q[1];
rz(2.9309978) q[3];
sx q[3];
rz(-0.67376332) q[3];
sx q[3];
rz(-1.3099561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.148968) q[2];
sx q[2];
rz(-0.31105369) q[2];
sx q[2];
rz(1.0827433) q[2];
rz(3.0454214) q[3];
sx q[3];
rz(-1.440719) q[3];
sx q[3];
rz(1.910803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(0.37725317) q[0];
sx q[0];
rz(-0.23582533) q[0];
sx q[0];
rz(2.8073231) q[0];
rz(-1.9175247) q[1];
sx q[1];
rz(-1.5609488) q[1];
sx q[1];
rz(2.8589378) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44336549) q[0];
sx q[0];
rz(-0.8015612) q[0];
sx q[0];
rz(-2.482224) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3547782) q[2];
sx q[2];
rz(-1.2470494) q[2];
sx q[2];
rz(0.87460364) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54284755) q[1];
sx q[1];
rz(-1.7464906) q[1];
sx q[1];
rz(0.41136841) q[1];
x q[2];
rz(1.71404) q[3];
sx q[3];
rz(-0.69056615) q[3];
sx q[3];
rz(2.7018202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9294372) q[2];
sx q[2];
rz(-1.9781457) q[2];
sx q[2];
rz(-2.4003417) q[2];
rz(0.50179982) q[3];
sx q[3];
rz(-1.8024249) q[3];
sx q[3];
rz(0.58399502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0989477) q[0];
sx q[0];
rz(-0.89530033) q[0];
sx q[0];
rz(-2.6328971) q[0];
rz(-0.11518654) q[1];
sx q[1];
rz(-2.4337264) q[1];
sx q[1];
rz(0.68181109) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17446974) q[0];
sx q[0];
rz(-2.5645064) q[0];
sx q[0];
rz(-1.2601) q[0];
rz(-pi) q[1];
rz(2.290906) q[2];
sx q[2];
rz(-0.88973532) q[2];
sx q[2];
rz(-1.016664) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9224285) q[1];
sx q[1];
rz(-2.104496) q[1];
sx q[1];
rz(2.1283172) q[1];
rz(-pi) q[2];
rz(0.82716771) q[3];
sx q[3];
rz(-0.70561545) q[3];
sx q[3];
rz(-0.08882113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.29601413) q[2];
sx q[2];
rz(-1.5073551) q[2];
sx q[2];
rz(-2.8005023) q[2];
rz(2.0579445) q[3];
sx q[3];
rz(-0.79569474) q[3];
sx q[3];
rz(-2.9705689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1682128) q[0];
sx q[0];
rz(-1.6642878) q[0];
sx q[0];
rz(-0.99933495) q[0];
rz(2.534261) q[1];
sx q[1];
rz(-2.2317531) q[1];
sx q[1];
rz(-2.9324525) q[1];
rz(2.4773459) q[2];
sx q[2];
rz(-1.9953809) q[2];
sx q[2];
rz(-2.8473163) q[2];
rz(2.7568983) q[3];
sx q[3];
rz(-1.7582498) q[3];
sx q[3];
rz(2.869217) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
