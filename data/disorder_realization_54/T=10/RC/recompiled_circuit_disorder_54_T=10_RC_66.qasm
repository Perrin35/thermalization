OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3611276) q[0];
sx q[0];
rz(-0.68929231) q[0];
sx q[0];
rz(-0.33049345) q[0];
rz(0.32980546) q[1];
sx q[1];
rz(5.4332241) q[1];
sx q[1];
rz(10.133893) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015895122) q[0];
sx q[0];
rz(-2.1791434) q[0];
sx q[0];
rz(-2.7817821) q[0];
rz(-pi) q[1];
rz(3.1094413) q[2];
sx q[2];
rz(-1.8841779) q[2];
sx q[2];
rz(-1.9762447) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3173649) q[1];
sx q[1];
rz(-2.1050276) q[1];
sx q[1];
rz(-0.51170106) q[1];
rz(-0.96592824) q[3];
sx q[3];
rz(-0.60797193) q[3];
sx q[3];
rz(1.0498429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7314529) q[2];
sx q[2];
rz(-2.7316766) q[2];
sx q[2];
rz(-1.6072134) q[2];
rz(2.2065227) q[3];
sx q[3];
rz(-1.3053223) q[3];
sx q[3];
rz(2.3527761) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4193029) q[0];
sx q[0];
rz(-3.0794444) q[0];
sx q[0];
rz(0.62227917) q[0];
rz(-0.17624804) q[1];
sx q[1];
rz(-1.9269678) q[1];
sx q[1];
rz(2.2252749) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30442552) q[0];
sx q[0];
rz(-0.93351782) q[0];
sx q[0];
rz(-0.62563719) q[0];
x q[1];
rz(2.3071438) q[2];
sx q[2];
rz(-1.4730028) q[2];
sx q[2];
rz(-2.9546839) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9333222) q[1];
sx q[1];
rz(-1.2332488) q[1];
sx q[1];
rz(-0.34668215) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5171492) q[3];
sx q[3];
rz(-0.39857736) q[3];
sx q[3];
rz(0.22565354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.24094412) q[2];
sx q[2];
rz(-0.99664656) q[2];
sx q[2];
rz(-0.82143482) q[2];
rz(-0.017283043) q[3];
sx q[3];
rz(-0.79764962) q[3];
sx q[3];
rz(-2.3582874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9538552) q[0];
sx q[0];
rz(-1.563235) q[0];
sx q[0];
rz(1.0082555) q[0];
rz(-3.1058274) q[1];
sx q[1];
rz(-1.6556031) q[1];
sx q[1];
rz(-2.6170513) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13428155) q[0];
sx q[0];
rz(-1.6005922) q[0];
sx q[0];
rz(1.7344463) q[0];
x q[1];
rz(1.8345941) q[2];
sx q[2];
rz(-2.2321777) q[2];
sx q[2];
rz(-1.7922572) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9959065) q[1];
sx q[1];
rz(-1.6325163) q[1];
sx q[1];
rz(-1.593959) q[1];
x q[2];
rz(-2.9933661) q[3];
sx q[3];
rz(-1.1550511) q[3];
sx q[3];
rz(1.1566597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.77164578) q[2];
sx q[2];
rz(-1.5180072) q[2];
sx q[2];
rz(-1.6463722) q[2];
rz(-1.3211936) q[3];
sx q[3];
rz(-1.4097872) q[3];
sx q[3];
rz(-2.7041919) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33048531) q[0];
sx q[0];
rz(-0.91902584) q[0];
sx q[0];
rz(3.0526429) q[0];
rz(-2.6308909) q[1];
sx q[1];
rz(-2.316244) q[1];
sx q[1];
rz(0.68960062) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8423858) q[0];
sx q[0];
rz(-0.33873522) q[0];
sx q[0];
rz(2.6329106) q[0];
rz(-pi) q[1];
rz(-1.3454516) q[2];
sx q[2];
rz(-1.0190522) q[2];
sx q[2];
rz(-1.6793041) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.00011132414) q[1];
sx q[1];
rz(-2.3173216) q[1];
sx q[1];
rz(-0.4517171) q[1];
rz(0.27130178) q[3];
sx q[3];
rz(-0.50008431) q[3];
sx q[3];
rz(0.015319583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4758063) q[2];
sx q[2];
rz(-2.4643354) q[2];
sx q[2];
rz(-0.62292567) q[2];
rz(-1.1359435) q[3];
sx q[3];
rz(-2.9562852) q[3];
sx q[3];
rz(-0.46666551) q[3];
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
rz(-2.2899807) q[0];
sx q[0];
rz(-1.8988134) q[0];
sx q[0];
rz(-0.583453) q[0];
rz(1.1460229) q[1];
sx q[1];
rz(-1.4910411) q[1];
sx q[1];
rz(1.6437644) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9153584) q[0];
sx q[0];
rz(-0.81904531) q[0];
sx q[0];
rz(0.35924964) q[0];
x q[1];
rz(1.6696817) q[2];
sx q[2];
rz(-2.3660198) q[2];
sx q[2];
rz(2.0688187) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.51753804) q[1];
sx q[1];
rz(-0.86746403) q[1];
sx q[1];
rz(-1.6066949) q[1];
rz(-pi) q[2];
rz(1.7329526) q[3];
sx q[3];
rz(-2.080653) q[3];
sx q[3];
rz(-1.7061403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4862165) q[2];
sx q[2];
rz(-1.5950173) q[2];
sx q[2];
rz(-1.5636469) q[2];
rz(-0.90562138) q[3];
sx q[3];
rz(-0.27035299) q[3];
sx q[3];
rz(-0.63846987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49801302) q[0];
sx q[0];
rz(-1.4587412) q[0];
sx q[0];
rz(-3.0946099) q[0];
rz(-2.9934096) q[1];
sx q[1];
rz(-2.2427185) q[1];
sx q[1];
rz(1.7061589) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2611321) q[0];
sx q[0];
rz(-1.4229703) q[0];
sx q[0];
rz(0.79974215) q[0];
rz(-pi) q[1];
rz(2.1204505) q[2];
sx q[2];
rz(-1.6209941) q[2];
sx q[2];
rz(0.93027885) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1286436) q[1];
sx q[1];
rz(-1.1638906) q[1];
sx q[1];
rz(3.0055771) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6834429) q[3];
sx q[3];
rz(-2.9023691) q[3];
sx q[3];
rz(0.84157543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0753714) q[2];
sx q[2];
rz(-1.0152638) q[2];
sx q[2];
rz(0.071468778) q[2];
rz(-1.6890769) q[3];
sx q[3];
rz(-1.444933) q[3];
sx q[3];
rz(2.215109) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8554095) q[0];
sx q[0];
rz(-1.3278642) q[0];
sx q[0];
rz(0.57762161) q[0];
rz(1.2795992) q[1];
sx q[1];
rz(-1.7956693) q[1];
sx q[1];
rz(1.0095899) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5492064) q[0];
sx q[0];
rz(-1.2388595) q[0];
sx q[0];
rz(0.50360002) q[0];
x q[1];
rz(-0.52092123) q[2];
sx q[2];
rz(-2.705057) q[2];
sx q[2];
rz(1.6233363) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2958888) q[1];
sx q[1];
rz(-1.0755952) q[1];
sx q[1];
rz(-1.1843029) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61543492) q[3];
sx q[3];
rz(-0.86343599) q[3];
sx q[3];
rz(-0.59188852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0916831) q[2];
sx q[2];
rz(-2.78806) q[2];
sx q[2];
rz(1.1996777) q[2];
rz(0.47232929) q[3];
sx q[3];
rz(-1.4940558) q[3];
sx q[3];
rz(-2.1388617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49884477) q[0];
sx q[0];
rz(-1.5002748) q[0];
sx q[0];
rz(-2.902466) q[0];
rz(0.7827951) q[1];
sx q[1];
rz(-2.6233964) q[1];
sx q[1];
rz(-2.696864) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.037576588) q[0];
sx q[0];
rz(-2.4009279) q[0];
sx q[0];
rz(0.91233493) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4773981) q[2];
sx q[2];
rz(-2.1014629) q[2];
sx q[2];
rz(-1.2150089) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4152894) q[1];
sx q[1];
rz(-0.88978926) q[1];
sx q[1];
rz(-2.9773832) q[1];
x q[2];
rz(-0.41110699) q[3];
sx q[3];
rz(-1.4484222) q[3];
sx q[3];
rz(-1.3947595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7198221) q[2];
sx q[2];
rz(-0.82227102) q[2];
sx q[2];
rz(-2.3262809) q[2];
rz(-2.1067965) q[3];
sx q[3];
rz(-1.5765604) q[3];
sx q[3];
rz(-1.8243272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3141044) q[0];
sx q[0];
rz(-2.1056471) q[0];
sx q[0];
rz(-0.28717336) q[0];
rz(2.9526967) q[1];
sx q[1];
rz(-0.72383988) q[1];
sx q[1];
rz(2.8093991) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9438292) q[0];
sx q[0];
rz(-1.257886) q[0];
sx q[0];
rz(0.80425941) q[0];
x q[1];
rz(1.5154778) q[2];
sx q[2];
rz(-2.4419867) q[2];
sx q[2];
rz(-0.31457065) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5933983) q[1];
sx q[1];
rz(-1.8264923) q[1];
sx q[1];
rz(-2.7138608) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6094749) q[3];
sx q[3];
rz(-1.9241153) q[3];
sx q[3];
rz(-1.6649099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2087848) q[2];
sx q[2];
rz(-2.1366426) q[2];
sx q[2];
rz(2.3020111) q[2];
rz(1.2906637) q[3];
sx q[3];
rz(-1.7920114) q[3];
sx q[3];
rz(-2.4850142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0861417) q[0];
sx q[0];
rz(-0.76776004) q[0];
sx q[0];
rz(1.0593876) q[0];
rz(-2.1620031) q[1];
sx q[1];
rz(-1.1971808) q[1];
sx q[1];
rz(-0.25451452) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018325018) q[0];
sx q[0];
rz(-1.3402481) q[0];
sx q[0];
rz(0.075320764) q[0];
rz(0.43912402) q[2];
sx q[2];
rz(-0.52999485) q[2];
sx q[2];
rz(2.8119171) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.49452457) q[1];
sx q[1];
rz(-1.2459323) q[1];
sx q[1];
rz(-1.4855794) q[1];
rz(-pi) q[2];
rz(2.9322846) q[3];
sx q[3];
rz(-2.6946687) q[3];
sx q[3];
rz(-0.79093864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0467726) q[2];
sx q[2];
rz(-2.5929055) q[2];
sx q[2];
rz(2.0937031) q[2];
rz(1.3607599) q[3];
sx q[3];
rz(-2.4436617) q[3];
sx q[3];
rz(-0.60539436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1142674) q[0];
sx q[0];
rz(-1.1189168) q[0];
sx q[0];
rz(3.0615321) q[0];
rz(2.7813773) q[1];
sx q[1];
rz(-1.4604912) q[1];
sx q[1];
rz(2.1866658) q[1];
rz(1.9706456) q[2];
sx q[2];
rz(-0.56213899) q[2];
sx q[2];
rz(3.0664372) q[2];
rz(1.249282) q[3];
sx q[3];
rz(-0.63890639) q[3];
sx q[3];
rz(-1.3756868) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
