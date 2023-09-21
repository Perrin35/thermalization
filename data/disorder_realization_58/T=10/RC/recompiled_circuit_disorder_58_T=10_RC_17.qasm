OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9392202) q[0];
sx q[0];
rz(-0.4063172) q[0];
sx q[0];
rz(-2.321474) q[0];
rz(-0.36110538) q[1];
sx q[1];
rz(-2.5087924) q[1];
sx q[1];
rz(-2.3109205) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3842073) q[0];
sx q[0];
rz(-1.877458) q[0];
sx q[0];
rz(2.7469809) q[0];
x q[1];
rz(-1.0933502) q[2];
sx q[2];
rz(-0.89077836) q[2];
sx q[2];
rz(0.29104656) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7100916) q[1];
sx q[1];
rz(-2.3570865) q[1];
sx q[1];
rz(2.0246519) q[1];
rz(-2.41483) q[3];
sx q[3];
rz(-0.57075497) q[3];
sx q[3];
rz(-1.4445514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.43710199) q[2];
sx q[2];
rz(-2.7097242) q[2];
sx q[2];
rz(0.1201771) q[2];
rz(-1.9834571) q[3];
sx q[3];
rz(-1.3985671) q[3];
sx q[3];
rz(2.2226298) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0607818) q[0];
sx q[0];
rz(-1.3588384) q[0];
sx q[0];
rz(-0.91180116) q[0];
rz(0.78951019) q[1];
sx q[1];
rz(-0.98840886) q[1];
sx q[1];
rz(2.8149014) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3347496) q[0];
sx q[0];
rz(-2.202889) q[0];
sx q[0];
rz(0.64081162) q[0];
rz(-pi) q[1];
rz(1.8961043) q[2];
sx q[2];
rz(-0.46519687) q[2];
sx q[2];
rz(2.543769) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8596526) q[1];
sx q[1];
rz(-1.9743866) q[1];
sx q[1];
rz(-2.2380026) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0558415) q[3];
sx q[3];
rz(-0.88480703) q[3];
sx q[3];
rz(0.12967295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6313173) q[2];
sx q[2];
rz(-2.3687506) q[2];
sx q[2];
rz(-1.2871683) q[2];
rz(-0.10989799) q[3];
sx q[3];
rz(-1.7307614) q[3];
sx q[3];
rz(1.3818285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54863769) q[0];
sx q[0];
rz(-2.4023963) q[0];
sx q[0];
rz(-2.8116995) q[0];
rz(-0.27711162) q[1];
sx q[1];
rz(-1.3169293) q[1];
sx q[1];
rz(-1.057391) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7085416) q[0];
sx q[0];
rz(-1.2305224) q[0];
sx q[0];
rz(0.021854594) q[0];
rz(1.7705275) q[2];
sx q[2];
rz(-1.9057416) q[2];
sx q[2];
rz(2.0470326) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8101465) q[1];
sx q[1];
rz(-2.3815037) q[1];
sx q[1];
rz(2.0600832) q[1];
x q[2];
rz(-1.7626761) q[3];
sx q[3];
rz(-1.9865611) q[3];
sx q[3];
rz(-0.92431812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3588336) q[2];
sx q[2];
rz(-1.1153355) q[2];
sx q[2];
rz(-1.7791629) q[2];
rz(2.5168915) q[3];
sx q[3];
rz(-1.0995068) q[3];
sx q[3];
rz(-0.81956285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7841004) q[0];
sx q[0];
rz(-0.52629137) q[0];
sx q[0];
rz(-2.6065361) q[0];
rz(-1.1401945) q[1];
sx q[1];
rz(-1.813872) q[1];
sx q[1];
rz(0.16539703) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3446942) q[0];
sx q[0];
rz(-1.7555153) q[0];
sx q[0];
rz(0.12634191) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9781614) q[2];
sx q[2];
rz(-1.7497517) q[2];
sx q[2];
rz(-2.2043259) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.0488102) q[1];
sx q[1];
rz(-0.65199344) q[1];
sx q[1];
rz(-0.55744967) q[1];
rz(-pi) q[2];
rz(-2.1256251) q[3];
sx q[3];
rz(-1.8053683) q[3];
sx q[3];
rz(-0.50333422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7455204) q[2];
sx q[2];
rz(-1.8053651) q[2];
sx q[2];
rz(-2.9361434) q[2];
rz(1.127634) q[3];
sx q[3];
rz(-1.1536359) q[3];
sx q[3];
rz(-1.2566465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66185343) q[0];
sx q[0];
rz(-0.91402188) q[0];
sx q[0];
rz(2.7868295) q[0];
rz(-1.154249) q[1];
sx q[1];
rz(-2.216979) q[1];
sx q[1];
rz(0.23194557) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35271586) q[0];
sx q[0];
rz(-2.0577601) q[0];
sx q[0];
rz(-2.6608174) q[0];
rz(-pi) q[1];
rz(1.4300214) q[2];
sx q[2];
rz(-0.88015926) q[2];
sx q[2];
rz(1.8292793) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9310589) q[1];
sx q[1];
rz(-1.3506883) q[1];
sx q[1];
rz(-1.7935497) q[1];
rz(-pi) q[2];
rz(0.050458126) q[3];
sx q[3];
rz(-2.7726463) q[3];
sx q[3];
rz(-1.3621393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.140124) q[2];
sx q[2];
rz(-1.8073558) q[2];
sx q[2];
rz(-1.9011964) q[2];
rz(0.59605789) q[3];
sx q[3];
rz(-1.3052992) q[3];
sx q[3];
rz(1.8381455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0444788) q[0];
sx q[0];
rz(-3.0713186) q[0];
sx q[0];
rz(0.20275673) q[0];
rz(0.98908201) q[1];
sx q[1];
rz(-1.6977856) q[1];
sx q[1];
rz(-2.1441377) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4076685) q[0];
sx q[0];
rz(-1.4364388) q[0];
sx q[0];
rz(-0.863048) q[0];
rz(1.9428271) q[2];
sx q[2];
rz(-0.87429201) q[2];
sx q[2];
rz(-0.55559413) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3012078) q[1];
sx q[1];
rz(-1.9921229) q[1];
sx q[1];
rz(2.8788484) q[1];
rz(-pi) q[2];
rz(2.0632083) q[3];
sx q[3];
rz(-0.78740722) q[3];
sx q[3];
rz(1.2315962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.22770195) q[2];
sx q[2];
rz(-1.9753549) q[2];
sx q[2];
rz(-2.690199) q[2];
rz(-0.40870062) q[3];
sx q[3];
rz(-1.5390076) q[3];
sx q[3];
rz(1.9394978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.30615) q[0];
sx q[0];
rz(-2.7375484) q[0];
sx q[0];
rz(2.5174482) q[0];
rz(-1.5165326) q[1];
sx q[1];
rz(-2.8846278) q[1];
sx q[1];
rz(2.5278032) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95425883) q[0];
sx q[0];
rz(-1.7566534) q[0];
sx q[0];
rz(-2.3128187) q[0];
x q[1];
rz(1.5067528) q[2];
sx q[2];
rz(-1.2876858) q[2];
sx q[2];
rz(1.7664906) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2144199) q[1];
sx q[1];
rz(-1.629429) q[1];
sx q[1];
rz(-0.10417948) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3404487) q[3];
sx q[3];
rz(-2.6655572) q[3];
sx q[3];
rz(1.7649094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7081786) q[2];
sx q[2];
rz(-2.2257979) q[2];
sx q[2];
rz(1.1748574) q[2];
rz(2.5332149) q[3];
sx q[3];
rz(-1.7374246) q[3];
sx q[3];
rz(1.7181989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2414918) q[0];
sx q[0];
rz(-1.8429723) q[0];
sx q[0];
rz(2.6532145) q[0];
rz(-1.6237367) q[1];
sx q[1];
rz(-1.7428215) q[1];
sx q[1];
rz(2.1571295) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9280633) q[0];
sx q[0];
rz(-2.5941879) q[0];
sx q[0];
rz(-0.071896032) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77163561) q[2];
sx q[2];
rz(-1.2727591) q[2];
sx q[2];
rz(-2.3723797) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.24482803) q[1];
sx q[1];
rz(-0.36537376) q[1];
sx q[1];
rz(0.016896292) q[1];
rz(-pi) q[2];
rz(0.55926178) q[3];
sx q[3];
rz(-0.54702938) q[3];
sx q[3];
rz(2.3032041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.70665923) q[2];
sx q[2];
rz(-1.6985396) q[2];
sx q[2];
rz(2.6064176) q[2];
rz(1.0501856) q[3];
sx q[3];
rz(-1.340056) q[3];
sx q[3];
rz(-1.0092658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.64131367) q[0];
sx q[0];
rz(-2.2135493) q[0];
sx q[0];
rz(-0.6859268) q[0];
rz(-0.39086875) q[1];
sx q[1];
rz(-2.1978244) q[1];
sx q[1];
rz(0.92591441) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5110916) q[0];
sx q[0];
rz(-1.4303659) q[0];
sx q[0];
rz(-1.7846084) q[0];
x q[1];
rz(-2.3472896) q[2];
sx q[2];
rz(-2.1199193) q[2];
sx q[2];
rz(1.6910451) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3013819) q[1];
sx q[1];
rz(-1.4242607) q[1];
sx q[1];
rz(-2.0135897) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9374398) q[3];
sx q[3];
rz(-1.6870058) q[3];
sx q[3];
rz(1.0128563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9459076) q[2];
sx q[2];
rz(-0.47161272) q[2];
sx q[2];
rz(-2.6055028) q[2];
rz(2.7219971) q[3];
sx q[3];
rz(-0.59046888) q[3];
sx q[3];
rz(1.6528116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-1.0855899) q[0];
sx q[0];
rz(-2.7828126) q[0];
sx q[0];
rz(-2.7465903) q[0];
rz(-1.6292054) q[1];
sx q[1];
rz(-1.869166) q[1];
sx q[1];
rz(-1.013247) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9625044) q[0];
sx q[0];
rz(-0.27217406) q[0];
sx q[0];
rz(-2.1129235) q[0];
rz(-0.35372325) q[2];
sx q[2];
rz(-1.612066) q[2];
sx q[2];
rz(1.8281787) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4229065) q[1];
sx q[1];
rz(-1.6675073) q[1];
sx q[1];
rz(-2.7194517) q[1];
rz(2.2447484) q[3];
sx q[3];
rz(-0.92771155) q[3];
sx q[3];
rz(-0.53900063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6932678) q[2];
sx q[2];
rz(-1.6222745) q[2];
sx q[2];
rz(0.75941336) q[2];
rz(1.365472) q[3];
sx q[3];
rz(-1.1080192) q[3];
sx q[3];
rz(0.56994462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7286745) q[0];
sx q[0];
rz(-1.8771794) q[0];
sx q[0];
rz(-1.9369453) q[0];
rz(0.57327523) q[1];
sx q[1];
rz(-1.9075867) q[1];
sx q[1];
rz(1.8214068) q[1];
rz(2.0157094) q[2];
sx q[2];
rz(-1.5013668) q[2];
sx q[2];
rz(-2.8253386) q[2];
rz(2.2686601) q[3];
sx q[3];
rz(-1.6193661) q[3];
sx q[3];
rz(-1.247874) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
