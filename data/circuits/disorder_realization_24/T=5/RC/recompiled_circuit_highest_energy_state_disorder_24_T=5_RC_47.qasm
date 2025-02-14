OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6744106) q[0];
sx q[0];
rz(-1.8888357) q[0];
sx q[0];
rz(2.3350265) q[0];
rz(-0.84713495) q[1];
sx q[1];
rz(-0.48589125) q[1];
sx q[1];
rz(-2.2810305) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0057324) q[0];
sx q[0];
rz(-1.4306698) q[0];
sx q[0];
rz(2.1707373) q[0];
x q[1];
rz(2.2769558) q[2];
sx q[2];
rz(-1.762391) q[2];
sx q[2];
rz(0.77048555) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.027801188) q[1];
sx q[1];
rz(-1.5209201) q[1];
sx q[1];
rz(-2.8667843) q[1];
rz(-0.99765649) q[3];
sx q[3];
rz(-1.1449779) q[3];
sx q[3];
rz(-0.65090685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.90031558) q[2];
sx q[2];
rz(-0.52854717) q[2];
sx q[2];
rz(-0.47145525) q[2];
rz(2.6491162) q[3];
sx q[3];
rz(-1.2329279) q[3];
sx q[3];
rz(-1.2537778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8660698) q[0];
sx q[0];
rz(-0.38551426) q[0];
sx q[0];
rz(-2.8634014) q[0];
rz(-1.1909852) q[1];
sx q[1];
rz(-1.8089801) q[1];
sx q[1];
rz(1.8461548) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4269276) q[0];
sx q[0];
rz(-0.97070995) q[0];
sx q[0];
rz(-0.99487181) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0619925) q[2];
sx q[2];
rz(-1.6222553) q[2];
sx q[2];
rz(2.1934794) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.55838096) q[1];
sx q[1];
rz(-1.19046) q[1];
sx q[1];
rz(-3.1391869) q[1];
rz(0.53923082) q[3];
sx q[3];
rz(-1.013621) q[3];
sx q[3];
rz(-2.2702366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1969455) q[2];
sx q[2];
rz(-1.6702007) q[2];
sx q[2];
rz(-2.6317281) q[2];
rz(-1.6950131) q[3];
sx q[3];
rz(-0.46049419) q[3];
sx q[3];
rz(-0.50486008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.533621) q[0];
sx q[0];
rz(-1.5190268) q[0];
sx q[0];
rz(0.98440379) q[0];
rz(0.016117485) q[1];
sx q[1];
rz(-1.1382444) q[1];
sx q[1];
rz(1.791753) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0371246) q[0];
sx q[0];
rz(-0.72513956) q[0];
sx q[0];
rz(2.8903236) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0764112) q[2];
sx q[2];
rz(-2.5177023) q[2];
sx q[2];
rz(2.9817493) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.74047725) q[1];
sx q[1];
rz(-0.7001895) q[1];
sx q[1];
rz(-3.0790331) q[1];
rz(2.3758395) q[3];
sx q[3];
rz(-1.4902876) q[3];
sx q[3];
rz(1.8395559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4237889) q[2];
sx q[2];
rz(-0.70222792) q[2];
sx q[2];
rz(0.88199893) q[2];
rz(-1.1841904) q[3];
sx q[3];
rz(-2.2004674) q[3];
sx q[3];
rz(0.73630303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0475567) q[0];
sx q[0];
rz(-1.2268257) q[0];
sx q[0];
rz(-3.0645698) q[0];
rz(-2.9432964) q[1];
sx q[1];
rz(-1.501333) q[1];
sx q[1];
rz(-2.95453) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3530981) q[0];
sx q[0];
rz(-1.9066296) q[0];
sx q[0];
rz(-1.737835) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18128975) q[2];
sx q[2];
rz(-1.224887) q[2];
sx q[2];
rz(-2.6225901) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.143972) q[1];
sx q[1];
rz(-1.5117497) q[1];
sx q[1];
rz(2.4322492) q[1];
rz(-pi) q[2];
rz(1.9105533) q[3];
sx q[3];
rz(-1.8700784) q[3];
sx q[3];
rz(1.2359985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5932811) q[2];
sx q[2];
rz(-2.686794) q[2];
sx q[2];
rz(1.6853257) q[2];
rz(0.71825394) q[3];
sx q[3];
rz(-1.3879489) q[3];
sx q[3];
rz(-0.83435241) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89183557) q[0];
sx q[0];
rz(-1.3351853) q[0];
sx q[0];
rz(0.43858132) q[0];
rz(-0.35722411) q[1];
sx q[1];
rz(-1.2780739) q[1];
sx q[1];
rz(-1.8908148) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8880576) q[0];
sx q[0];
rz(-1.6707289) q[0];
sx q[0];
rz(-1.9991831) q[0];
rz(-pi) q[1];
x q[1];
rz(0.078089291) q[2];
sx q[2];
rz(-1.2603344) q[2];
sx q[2];
rz(3.1333528) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.67235302) q[1];
sx q[1];
rz(-1.9682944) q[1];
sx q[1];
rz(2.7356477) q[1];
x q[2];
rz(-0.46535551) q[3];
sx q[3];
rz(-2.8520003) q[3];
sx q[3];
rz(-0.48894879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.232051) q[2];
sx q[2];
rz(-2.6614058) q[2];
sx q[2];
rz(1.5117744) q[2];
rz(3.1253641) q[3];
sx q[3];
rz(-2.0461693) q[3];
sx q[3];
rz(-2.3487263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69415724) q[0];
sx q[0];
rz(-1.8614391) q[0];
sx q[0];
rz(-1.9919027) q[0];
rz(-2.7206874) q[1];
sx q[1];
rz(-1.6890539) q[1];
sx q[1];
rz(3.0973869) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71647787) q[0];
sx q[0];
rz(-1.8938066) q[0];
sx q[0];
rz(-0.72854211) q[0];
x q[1];
rz(0.8755336) q[2];
sx q[2];
rz(-0.6630156) q[2];
sx q[2];
rz(-0.11817486) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.845712) q[1];
sx q[1];
rz(-2.7428877) q[1];
sx q[1];
rz(-1.7804342) q[1];
rz(-0.44852528) q[3];
sx q[3];
rz(-2.5618346) q[3];
sx q[3];
rz(3.1223402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9629024) q[2];
sx q[2];
rz(-0.89683214) q[2];
sx q[2];
rz(-2.1964591) q[2];
rz(-1.6804228) q[3];
sx q[3];
rz(-1.9943359) q[3];
sx q[3];
rz(-2.6233961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1951676) q[0];
sx q[0];
rz(-2.8312046) q[0];
sx q[0];
rz(2.2921966) q[0];
rz(1.7646344) q[1];
sx q[1];
rz(-0.57146776) q[1];
sx q[1];
rz(1.2845385) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021668108) q[0];
sx q[0];
rz(-2.4942022) q[0];
sx q[0];
rz(0.68418087) q[0];
rz(0.30751245) q[2];
sx q[2];
rz(-1.4685681) q[2];
sx q[2];
rz(-0.44164613) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0459874) q[1];
sx q[1];
rz(-0.52667499) q[1];
sx q[1];
rz(3.0707568) q[1];
x q[2];
rz(0.86817812) q[3];
sx q[3];
rz(-0.68553151) q[3];
sx q[3];
rz(0.037151046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.66594243) q[2];
sx q[2];
rz(-1.7721756) q[2];
sx q[2];
rz(-1.5671889) q[2];
rz(-0.33268467) q[3];
sx q[3];
rz(-1.0654457) q[3];
sx q[3];
rz(-2.1596597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5653) q[0];
sx q[0];
rz(-1.1197634) q[0];
sx q[0];
rz(2.0530307) q[0];
rz(0.92900485) q[1];
sx q[1];
rz(-1.6477511) q[1];
sx q[1];
rz(-2.8378024) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32919183) q[0];
sx q[0];
rz(-0.69965967) q[0];
sx q[0];
rz(-0.60075642) q[0];
x q[1];
rz(-2.1985198) q[2];
sx q[2];
rz(-0.20096261) q[2];
sx q[2];
rz(-2.497884) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9907797) q[1];
sx q[1];
rz(-2.7013198) q[1];
sx q[1];
rz(1.7899465) q[1];
rz(-pi) q[2];
rz(-2.5098544) q[3];
sx q[3];
rz(-2.4186385) q[3];
sx q[3];
rz(2.637104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2527689) q[2];
sx q[2];
rz(-1.4759651) q[2];
sx q[2];
rz(-2.8775173) q[2];
rz(1.1325599) q[3];
sx q[3];
rz(-2.8045636) q[3];
sx q[3];
rz(-2.0866709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7842512) q[0];
sx q[0];
rz(-0.34264523) q[0];
sx q[0];
rz(1.0145048) q[0];
rz(-2.2134589) q[1];
sx q[1];
rz(-1.4200297) q[1];
sx q[1];
rz(-2.6511505) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7785085) q[0];
sx q[0];
rz(-0.16101232) q[0];
sx q[0];
rz(1.9912316) q[0];
rz(-pi) q[1];
rz(1.1855256) q[2];
sx q[2];
rz(-1.8920355) q[2];
sx q[2];
rz(2.3209077) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0874493) q[1];
sx q[1];
rz(-1.6119401) q[1];
sx q[1];
rz(-1.270134) q[1];
rz(-pi) q[2];
rz(-2.1098299) q[3];
sx q[3];
rz(-1.7787053) q[3];
sx q[3];
rz(-2.5202005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8162615) q[2];
sx q[2];
rz(-1.0518495) q[2];
sx q[2];
rz(3.0391147) q[2];
rz(-1.8898194) q[3];
sx q[3];
rz(-1.7779558) q[3];
sx q[3];
rz(1.6583091) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3935299) q[0];
sx q[0];
rz(-0.04638014) q[0];
sx q[0];
rz(-1.8512132) q[0];
rz(-2.6875467) q[1];
sx q[1];
rz(-1.8171277) q[1];
sx q[1];
rz(-2.9581199) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1114706) q[0];
sx q[0];
rz(-2.9577177) q[0];
sx q[0];
rz(-1.4669815) q[0];
rz(-pi) q[1];
rz(3.119629) q[2];
sx q[2];
rz(-0.59810592) q[2];
sx q[2];
rz(-0.22207161) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.80663348) q[1];
sx q[1];
rz(-2.0628961) q[1];
sx q[1];
rz(-1.5978807) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0504543) q[3];
sx q[3];
rz(-2.1347858) q[3];
sx q[3];
rz(0.48479776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4219249) q[2];
sx q[2];
rz(-1.0591966) q[2];
sx q[2];
rz(0.025040778) q[2];
rz(2.5981564) q[3];
sx q[3];
rz(-0.70356026) q[3];
sx q[3];
rz(-2.8652625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80171361) q[0];
sx q[0];
rz(-2.1262953) q[0];
sx q[0];
rz(-2.3332818) q[0];
rz(0.33357757) q[1];
sx q[1];
rz(-1.7671276) q[1];
sx q[1];
rz(-0.50552013) q[1];
rz(0.12814604) q[2];
sx q[2];
rz(-1.2432877) q[2];
sx q[2];
rz(0.23040085) q[2];
rz(0.65350914) q[3];
sx q[3];
rz(-2.1353888) q[3];
sx q[3];
rz(0.63513811) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
