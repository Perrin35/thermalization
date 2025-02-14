OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0754492) q[0];
sx q[0];
rz(-1.0439405) q[0];
sx q[0];
rz(-3.1312842) q[0];
rz(2.161624) q[1];
sx q[1];
rz(-1.4449395) q[1];
sx q[1];
rz(-0.57656062) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2808997) q[0];
sx q[0];
rz(-1.4497978) q[0];
sx q[0];
rz(2.8314986) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61277436) q[2];
sx q[2];
rz(-0.98721993) q[2];
sx q[2];
rz(0.60223168) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.7652055) q[1];
sx q[1];
rz(-0.55843267) q[1];
sx q[1];
rz(-2.8064578) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.16657942) q[3];
sx q[3];
rz(-0.84462476) q[3];
sx q[3];
rz(2.8960901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6171241) q[2];
sx q[2];
rz(-1.4602129) q[2];
sx q[2];
rz(-1.8560393) q[2];
rz(1.5517976) q[3];
sx q[3];
rz(-1.0714622) q[3];
sx q[3];
rz(-2.0901399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1262421) q[0];
sx q[0];
rz(-1.9168357) q[0];
sx q[0];
rz(2.8958564) q[0];
rz(-1.0579146) q[1];
sx q[1];
rz(-1.825288) q[1];
sx q[1];
rz(2.8071075) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6049371) q[0];
sx q[0];
rz(-0.70022445) q[0];
sx q[0];
rz(0.74857736) q[0];
rz(-pi) q[1];
rz(1.9742786) q[2];
sx q[2];
rz(-2.2902787) q[2];
sx q[2];
rz(-0.011842273) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1002754) q[1];
sx q[1];
rz(-1.9875437) q[1];
sx q[1];
rz(-0.4112501) q[1];
rz(2.7354419) q[3];
sx q[3];
rz(-2.5037919) q[3];
sx q[3];
rz(0.33406182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1841396) q[2];
sx q[2];
rz(-0.89749557) q[2];
sx q[2];
rz(1.755836) q[2];
rz(-2.1615248) q[3];
sx q[3];
rz(-0.46531427) q[3];
sx q[3];
rz(0.050962713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7053213) q[0];
sx q[0];
rz(-0.956981) q[0];
sx q[0];
rz(1.7514239) q[0];
rz(0.35274371) q[1];
sx q[1];
rz(-2.2309062) q[1];
sx q[1];
rz(0.62044755) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.432029) q[0];
sx q[0];
rz(-1.480446) q[0];
sx q[0];
rz(-1.4960775) q[0];
rz(-pi) q[1];
rz(-2.1587055) q[2];
sx q[2];
rz(-1.6574142) q[2];
sx q[2];
rz(-0.69332214) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7812278) q[1];
sx q[1];
rz(-1.2457799) q[1];
sx q[1];
rz(-2.8311391) q[1];
x q[2];
rz(-2.6246715) q[3];
sx q[3];
rz(-2.1955588) q[3];
sx q[3];
rz(-2.3464835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0188401) q[2];
sx q[2];
rz(-2.7285125) q[2];
sx q[2];
rz(1.2472461) q[2];
rz(2.9774169) q[3];
sx q[3];
rz(-1.7069867) q[3];
sx q[3];
rz(-2.5119787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.71488798) q[0];
sx q[0];
rz(-0.96788228) q[0];
sx q[0];
rz(-1.2870652) q[0];
rz(0.60028589) q[1];
sx q[1];
rz(-1.7900107) q[1];
sx q[1];
rz(0.90369019) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1335588) q[0];
sx q[0];
rz(-1.6370336) q[0];
sx q[0];
rz(-3.0341107) q[0];
rz(-pi) q[1];
rz(2.7616384) q[2];
sx q[2];
rz(-1.5075052) q[2];
sx q[2];
rz(1.5587057) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7052482) q[1];
sx q[1];
rz(-1.3911934) q[1];
sx q[1];
rz(-0.257538) q[1];
rz(-pi) q[2];
rz(-2.4900808) q[3];
sx q[3];
rz(-2.6389696) q[3];
sx q[3];
rz(-2.9732957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51957447) q[2];
sx q[2];
rz(-0.833424) q[2];
sx q[2];
rz(2.3681417) q[2];
rz(1.8526239) q[3];
sx q[3];
rz(-0.52923146) q[3];
sx q[3];
rz(2.4066063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89161038) q[0];
sx q[0];
rz(-2.1365428) q[0];
sx q[0];
rz(-2.6494359) q[0];
rz(2.6026169) q[1];
sx q[1];
rz(-1.69918) q[1];
sx q[1];
rz(2.0416562) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5163706) q[0];
sx q[0];
rz(-1.8305792) q[0];
sx q[0];
rz(-2.7894839) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1824838) q[2];
sx q[2];
rz(-1.5619742) q[2];
sx q[2];
rz(-0.1977284) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.392994) q[1];
sx q[1];
rz(-1.4778607) q[1];
sx q[1];
rz(0.8108906) q[1];
rz(-pi) q[2];
rz(1.8995011) q[3];
sx q[3];
rz(-2.0315758) q[3];
sx q[3];
rz(-0.24390175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.34053549) q[2];
sx q[2];
rz(-0.33581442) q[2];
sx q[2];
rz(-0.56279969) q[2];
rz(2.4417012) q[3];
sx q[3];
rz(-1.8507345) q[3];
sx q[3];
rz(-2.8904397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7971147) q[0];
sx q[0];
rz(-0.7239224) q[0];
sx q[0];
rz(0.35514721) q[0];
rz(-1.2190602) q[1];
sx q[1];
rz(-1.2390169) q[1];
sx q[1];
rz(-0.452279) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41846965) q[0];
sx q[0];
rz(-2.4112066) q[0];
sx q[0];
rz(2.1900221) q[0];
rz(-pi) q[1];
rz(0.81461774) q[2];
sx q[2];
rz(-2.6067408) q[2];
sx q[2];
rz(1.2592821) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2792795) q[1];
sx q[1];
rz(-1.24839) q[1];
sx q[1];
rz(0.025397852) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5127453) q[3];
sx q[3];
rz(-1.4553242) q[3];
sx q[3];
rz(-2.2001298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5111115) q[2];
sx q[2];
rz(-2.6177572) q[2];
sx q[2];
rz(-3.0044978) q[2];
rz(-1.5824205) q[3];
sx q[3];
rz(-1.987792) q[3];
sx q[3];
rz(-0.77663511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95098507) q[0];
sx q[0];
rz(-0.76222104) q[0];
sx q[0];
rz(-1.5451587) q[0];
rz(0.51721382) q[1];
sx q[1];
rz(-1.9904174) q[1];
sx q[1];
rz(2.1545765) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0407437) q[0];
sx q[0];
rz(-3.0090539) q[0];
sx q[0];
rz(-2.9805471) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1783243) q[2];
sx q[2];
rz(-1.844256) q[2];
sx q[2];
rz(0.69404049) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4728022) q[1];
sx q[1];
rz(-1.3061532) q[1];
sx q[1];
rz(1.8084779) q[1];
x q[2];
rz(-1.2521947) q[3];
sx q[3];
rz(-1.650839) q[3];
sx q[3];
rz(0.53832507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3551657) q[2];
sx q[2];
rz(-2.1106909) q[2];
sx q[2];
rz(-0.008358566) q[2];
rz(-2.9110294) q[3];
sx q[3];
rz(-1.2059261) q[3];
sx q[3];
rz(1.4768538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.01345988) q[0];
sx q[0];
rz(-3.1210493) q[0];
sx q[0];
rz(1.531456) q[0];
rz(1.5229335) q[1];
sx q[1];
rz(-1.5956655) q[1];
sx q[1];
rz(-1.5628372) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7195936) q[0];
sx q[0];
rz(-2.5223456) q[0];
sx q[0];
rz(-0.25281711) q[0];
rz(-pi) q[1];
rz(-1.426307) q[2];
sx q[2];
rz(-1.943271) q[2];
sx q[2];
rz(2.5691751) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.48846515) q[1];
sx q[1];
rz(-1.9886977) q[1];
sx q[1];
rz(-2.6326219) q[1];
rz(-pi) q[2];
rz(2.7430629) q[3];
sx q[3];
rz(-0.70793286) q[3];
sx q[3];
rz(-1.7205451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.83850399) q[2];
sx q[2];
rz(-1.9946626) q[2];
sx q[2];
rz(3.1040891) q[2];
rz(0.23669067) q[3];
sx q[3];
rz(-2.762837) q[3];
sx q[3];
rz(1.8481351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26838747) q[0];
sx q[0];
rz(-2.3516042) q[0];
sx q[0];
rz(2.068212) q[0];
rz(2.7997596) q[1];
sx q[1];
rz(-2.5624202) q[1];
sx q[1];
rz(2.5198708) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8499924) q[0];
sx q[0];
rz(-1.1860285) q[0];
sx q[0];
rz(-1.8509393) q[0];
rz(-1.4474117) q[2];
sx q[2];
rz(-0.91431352) q[2];
sx q[2];
rz(-0.8379762) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.016704233) q[1];
sx q[1];
rz(-2.0604265) q[1];
sx q[1];
rz(-1.917065) q[1];
rz(-pi) q[2];
rz(0.55523606) q[3];
sx q[3];
rz(-0.16517565) q[3];
sx q[3];
rz(0.092008807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6217893) q[2];
sx q[2];
rz(-2.9534464) q[2];
sx q[2];
rz(0.56358799) q[2];
rz(1.8300736) q[3];
sx q[3];
rz(-1.2389641) q[3];
sx q[3];
rz(0.81609503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9594864) q[0];
sx q[0];
rz(-0.86903787) q[0];
sx q[0];
rz(2.2667789) q[0];
rz(-1.4080217) q[1];
sx q[1];
rz(-2.4751016) q[1];
sx q[1];
rz(-0.63546884) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59171133) q[0];
sx q[0];
rz(-2.3259865) q[0];
sx q[0];
rz(2.7278103) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9866139) q[2];
sx q[2];
rz(-1.8228616) q[2];
sx q[2];
rz(-0.76463715) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.26837743) q[1];
sx q[1];
rz(-2.4627373) q[1];
sx q[1];
rz(-2.2427903) q[1];
x q[2];
rz(-2.1346774) q[3];
sx q[3];
rz(-1.0936519) q[3];
sx q[3];
rz(-1.7686219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.94329876) q[2];
sx q[2];
rz(-1.0903) q[2];
sx q[2];
rz(2.3401006) q[2];
rz(-1.21579) q[3];
sx q[3];
rz(-0.91167584) q[3];
sx q[3];
rz(3.099814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53032482) q[0];
sx q[0];
rz(-1.7014736) q[0];
sx q[0];
rz(0.3076719) q[0];
rz(-1.2904185) q[1];
sx q[1];
rz(-2.1967874) q[1];
sx q[1];
rz(-2.8312942) q[1];
rz(-0.18193131) q[2];
sx q[2];
rz(-1.4375189) q[2];
sx q[2];
rz(0.8663485) q[2];
rz(-2.484816) q[3];
sx q[3];
rz(-1.4639502) q[3];
sx q[3];
rz(-1.7456036) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
