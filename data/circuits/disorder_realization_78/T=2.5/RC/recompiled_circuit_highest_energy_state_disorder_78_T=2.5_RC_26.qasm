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
rz(-3.0566293) q[0];
sx q[0];
rz(-0.30240107) q[0];
sx q[0];
rz(0.032057134) q[0];
rz(-0.0078460296) q[1];
sx q[1];
rz(-0.50385952) q[1];
sx q[1];
rz(-0.75262466) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0085474) q[0];
sx q[0];
rz(-0.38131443) q[0];
sx q[0];
rz(-0.56579907) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7787894) q[2];
sx q[2];
rz(-1.2191091) q[2];
sx q[2];
rz(1.9039409) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6258273) q[1];
sx q[1];
rz(-1.9099565) q[1];
sx q[1];
rz(1.304342) q[1];
rz(0.21423046) q[3];
sx q[3];
rz(-1.2204683) q[3];
sx q[3];
rz(-0.080834724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.96693119) q[2];
sx q[2];
rz(-1.9661247) q[2];
sx q[2];
rz(0.82007972) q[2];
rz(0.44102937) q[3];
sx q[3];
rz(-2.0377772) q[3];
sx q[3];
rz(2.5681514) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012861982) q[0];
sx q[0];
rz(-2.2282889) q[0];
sx q[0];
rz(-0.41020694) q[0];
rz(2.7606616) q[1];
sx q[1];
rz(-1.610264) q[1];
sx q[1];
rz(-1.358323) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5986431) q[0];
sx q[0];
rz(-0.99587593) q[0];
sx q[0];
rz(-0.47358124) q[0];
rz(0.50308268) q[2];
sx q[2];
rz(-1.0840992) q[2];
sx q[2];
rz(-2.8622735) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0002901) q[1];
sx q[1];
rz(-1.7385671) q[1];
sx q[1];
rz(3.0649158) q[1];
rz(-pi) q[2];
rz(-1.7776599) q[3];
sx q[3];
rz(-2.4914722) q[3];
sx q[3];
rz(0.93322414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7831948) q[2];
sx q[2];
rz(-0.40893778) q[2];
sx q[2];
rz(-0.50296339) q[2];
rz(-1.8424235) q[3];
sx q[3];
rz(-2.3830569) q[3];
sx q[3];
rz(2.9097596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7314887) q[0];
sx q[0];
rz(-0.51173156) q[0];
sx q[0];
rz(2.1657535) q[0];
rz(2.2052235) q[1];
sx q[1];
rz(-1.5872637) q[1];
sx q[1];
rz(0.13793129) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5463211) q[0];
sx q[0];
rz(-2.0450651) q[0];
sx q[0];
rz(2.7300937) q[0];
x q[1];
rz(0.39769002) q[2];
sx q[2];
rz(-1.4639336) q[2];
sx q[2];
rz(-1.7025089) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5974849) q[1];
sx q[1];
rz(-0.710809) q[1];
sx q[1];
rz(1.7729458) q[1];
rz(-pi) q[2];
rz(2.6343752) q[3];
sx q[3];
rz(-2.6196369) q[3];
sx q[3];
rz(1.0591266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5522573) q[2];
sx q[2];
rz(-1.6103123) q[2];
sx q[2];
rz(1.4790347) q[2];
rz(-0.79834437) q[3];
sx q[3];
rz(-0.84404293) q[3];
sx q[3];
rz(0.020309694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.289157) q[0];
sx q[0];
rz(-0.11141369) q[0];
sx q[0];
rz(-2.074746) q[0];
rz(1.5929219) q[1];
sx q[1];
rz(-1.8916062) q[1];
sx q[1];
rz(0.41935316) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33808483) q[0];
sx q[0];
rz(-0.28059059) q[0];
sx q[0];
rz(-2.1795033) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6134562) q[2];
sx q[2];
rz(-1.6491051) q[2];
sx q[2];
rz(2.9700043) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6450582) q[1];
sx q[1];
rz(-2.7784111) q[1];
sx q[1];
rz(-1.3892322) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0556424) q[3];
sx q[3];
rz(-1.8797698) q[3];
sx q[3];
rz(-0.072871836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7833917) q[2];
sx q[2];
rz(-2.6322067) q[2];
sx q[2];
rz(-1.5979213) q[2];
rz(1.915043) q[3];
sx q[3];
rz(-1.2164601) q[3];
sx q[3];
rz(-1.8515057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.876038) q[0];
sx q[0];
rz(-2.4198678) q[0];
sx q[0];
rz(0.75620404) q[0];
rz(-1.8887695) q[1];
sx q[1];
rz(-0.93106657) q[1];
sx q[1];
rz(1.4160215) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6373581) q[0];
sx q[0];
rz(-1.7880482) q[0];
sx q[0];
rz(0.603032) q[0];
rz(-pi) q[1];
rz(0.3438832) q[2];
sx q[2];
rz(-1.4522219) q[2];
sx q[2];
rz(-0.32394513) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.51235837) q[1];
sx q[1];
rz(-0.53881379) q[1];
sx q[1];
rz(0.013756559) q[1];
rz(-pi) q[2];
rz(2.6719499) q[3];
sx q[3];
rz(-2.7386463) q[3];
sx q[3];
rz(-2.6279891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8004134) q[2];
sx q[2];
rz(-0.74113733) q[2];
sx q[2];
rz(-0.18079147) q[2];
rz(1.4888658) q[3];
sx q[3];
rz(-1.3988262) q[3];
sx q[3];
rz(1.5159877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7086696) q[0];
sx q[0];
rz(-2.0661418) q[0];
sx q[0];
rz(2.3413626) q[0];
rz(2.8569787) q[1];
sx q[1];
rz(-1.0601284) q[1];
sx q[1];
rz(-3.0175993) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5381966) q[0];
sx q[0];
rz(-1.5423703) q[0];
sx q[0];
rz(-1.4409426) q[0];
rz(-pi) q[1];
rz(-1.8058067) q[2];
sx q[2];
rz(-1.5127276) q[2];
sx q[2];
rz(-0.11561671) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4912419) q[1];
sx q[1];
rz(-1.6524775) q[1];
sx q[1];
rz(0.30568576) q[1];
rz(-pi) q[2];
rz(1.049253) q[3];
sx q[3];
rz(-2.4774356) q[3];
sx q[3];
rz(-0.60598999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.93512145) q[2];
sx q[2];
rz(-1.4149041) q[2];
sx q[2];
rz(-0.07621152) q[2];
rz(1.8501836) q[3];
sx q[3];
rz(-1.094787) q[3];
sx q[3];
rz(-2.5230303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16977075) q[0];
sx q[0];
rz(-1.562028) q[0];
sx q[0];
rz(0.75702697) q[0];
rz(-0.79611671) q[1];
sx q[1];
rz(-1.7346104) q[1];
sx q[1];
rz(0.98181358) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9011089) q[0];
sx q[0];
rz(-2.3496858) q[0];
sx q[0];
rz(-1.9740941) q[0];
x q[1];
rz(-2.1481291) q[2];
sx q[2];
rz(-3.0138123) q[2];
sx q[2];
rz(1.8005231) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0701778) q[1];
sx q[1];
rz(-1.1648253) q[1];
sx q[1];
rz(-1.0557879) q[1];
rz(-pi) q[2];
rz(-1.8514093) q[3];
sx q[3];
rz(-1.7095437) q[3];
sx q[3];
rz(2.8900103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7299399) q[2];
sx q[2];
rz(-0.81092683) q[2];
sx q[2];
rz(0.5438424) q[2];
rz(-3.1206711) q[3];
sx q[3];
rz(-1.7048416) q[3];
sx q[3];
rz(0.81834832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92886096) q[0];
sx q[0];
rz(-0.91101557) q[0];
sx q[0];
rz(-0.62335706) q[0];
rz(2.8202672) q[1];
sx q[1];
rz(-2.5265381) q[1];
sx q[1];
rz(-0.26434937) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67548114) q[0];
sx q[0];
rz(-2.3947885) q[0];
sx q[0];
rz(1.8854499) q[0];
x q[1];
rz(-1.2351722) q[2];
sx q[2];
rz(-1.610029) q[2];
sx q[2];
rz(-0.36612636) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0496638) q[1];
sx q[1];
rz(-0.85159066) q[1];
sx q[1];
rz(2.767241) q[1];
x q[2];
rz(1.8885625) q[3];
sx q[3];
rz(-1.4444286) q[3];
sx q[3];
rz(1.7333836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.136772) q[2];
sx q[2];
rz(-2.4592168) q[2];
sx q[2];
rz(-0.47478673) q[2];
rz(2.7759806) q[3];
sx q[3];
rz(-0.48706278) q[3];
sx q[3];
rz(-0.18317187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8356165) q[0];
sx q[0];
rz(-2.7662179) q[0];
sx q[0];
rz(-2.6557652) q[0];
rz(-1.0659418) q[1];
sx q[1];
rz(-1.7168047) q[1];
sx q[1];
rz(0.33445439) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36223672) q[0];
sx q[0];
rz(-1.4051172) q[0];
sx q[0];
rz(-0.30596531) q[0];
rz(-pi) q[1];
rz(-2.8489001) q[2];
sx q[2];
rz(-0.45537696) q[2];
sx q[2];
rz(0.51560452) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0171256) q[1];
sx q[1];
rz(-1.4216058) q[1];
sx q[1];
rz(-0.45254947) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79631898) q[3];
sx q[3];
rz(-1.9069172) q[3];
sx q[3];
rz(1.6567865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7740384) q[2];
sx q[2];
rz(-2.2944821) q[2];
sx q[2];
rz(-0.96662194) q[2];
rz(-1.4878081) q[3];
sx q[3];
rz(-0.32854587) q[3];
sx q[3];
rz(-0.070913471) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8129355) q[0];
sx q[0];
rz(-0.85449496) q[0];
sx q[0];
rz(-0.27035126) q[0];
rz(1.6290172) q[1];
sx q[1];
rz(-0.83643475) q[1];
sx q[1];
rz(-0.62634748) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.644076) q[0];
sx q[0];
rz(-0.77267161) q[0];
sx q[0];
rz(2.0324043) q[0];
x q[1];
rz(-0.93339351) q[2];
sx q[2];
rz(-1.9248171) q[2];
sx q[2];
rz(0.4590946) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4463908) q[1];
sx q[1];
rz(-3.0425515) q[1];
sx q[1];
rz(0.42548577) q[1];
x q[2];
rz(2.7553431) q[3];
sx q[3];
rz(-1.874141) q[3];
sx q[3];
rz(1.4331499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7865929) q[2];
sx q[2];
rz(-2.2248416) q[2];
sx q[2];
rz(-1.5691441) q[2];
rz(2.3908424) q[3];
sx q[3];
rz(-0.34199491) q[3];
sx q[3];
rz(0.87740889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56562051) q[0];
sx q[0];
rz(-1.8436057) q[0];
sx q[0];
rz(2.5634503) q[0];
rz(1.7987953) q[1];
sx q[1];
rz(-0.31675757) q[1];
sx q[1];
rz(-2.996179) q[1];
rz(-3.0192791) q[2];
sx q[2];
rz(-1.5707471) q[2];
sx q[2];
rz(-1.6483501) q[2];
rz(0.20082898) q[3];
sx q[3];
rz(-2.8786537) q[3];
sx q[3];
rz(-2.636551) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
