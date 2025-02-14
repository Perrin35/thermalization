OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.00081113022) q[0];
sx q[0];
rz(-0.82511628) q[0];
sx q[0];
rz(0.42186475) q[0];
rz(-0.74692625) q[1];
sx q[1];
rz(3.7318228) q[1];
sx q[1];
rz(11.745315) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79347389) q[0];
sx q[0];
rz(-3.0464026) q[0];
sx q[0];
rz(-2.6497255) q[0];
rz(-pi) q[1];
rz(-2.9157588) q[2];
sx q[2];
rz(-1.3703114) q[2];
sx q[2];
rz(0.12910138) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5297048) q[1];
sx q[1];
rz(-0.48445736) q[1];
sx q[1];
rz(-0.92682181) q[1];
x q[2];
rz(1.2263359) q[3];
sx q[3];
rz(-1.4446961) q[3];
sx q[3];
rz(-0.9094224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.39516285) q[2];
sx q[2];
rz(-0.91150993) q[2];
sx q[2];
rz(2.5968623) q[2];
rz(1.644545) q[3];
sx q[3];
rz(-2.0293197) q[3];
sx q[3];
rz(1.3397608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6853952) q[0];
sx q[0];
rz(-2.39769) q[0];
sx q[0];
rz(0.5734545) q[0];
rz(0.042512976) q[1];
sx q[1];
rz(-1.3923693) q[1];
sx q[1];
rz(-1.2826077) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79270482) q[0];
sx q[0];
rz(-3.0137339) q[0];
sx q[0];
rz(2.9793745) q[0];
rz(-pi) q[1];
rz(1.3147215) q[2];
sx q[2];
rz(-1.8559858) q[2];
sx q[2];
rz(1.8023173) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8136422) q[1];
sx q[1];
rz(-0.54681603) q[1];
sx q[1];
rz(-0.27346404) q[1];
rz(-pi) q[2];
rz(-2.1759791) q[3];
sx q[3];
rz(-2.11016) q[3];
sx q[3];
rz(2.2210647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.87806988) q[2];
sx q[2];
rz(-1.9922549) q[2];
sx q[2];
rz(-2.6980706) q[2];
rz(-1.5524607) q[3];
sx q[3];
rz(-2.7510721) q[3];
sx q[3];
rz(-2.922557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012101128) q[0];
sx q[0];
rz(-0.87915593) q[0];
sx q[0];
rz(-0.40170676) q[0];
rz(1.7237639) q[1];
sx q[1];
rz(-2.6938853) q[1];
sx q[1];
rz(-1.1161425) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8153977) q[0];
sx q[0];
rz(-0.92741637) q[0];
sx q[0];
rz(-1.1872435) q[0];
rz(-pi) q[1];
rz(-0.1537519) q[2];
sx q[2];
rz(-0.88482403) q[2];
sx q[2];
rz(-1.5123715) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9520719) q[1];
sx q[1];
rz(-2.4671989) q[1];
sx q[1];
rz(0.86878573) q[1];
rz(-1.6468372) q[3];
sx q[3];
rz(-1.0234939) q[3];
sx q[3];
rz(-0.28441498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9125354) q[2];
sx q[2];
rz(-2.8450862) q[2];
sx q[2];
rz(-1.0079481) q[2];
rz(-1.6352765) q[3];
sx q[3];
rz(-1.9624036) q[3];
sx q[3];
rz(-2.262871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85663831) q[0];
sx q[0];
rz(-0.027218787) q[0];
sx q[0];
rz(-0.83754367) q[0];
rz(1.8800927) q[1];
sx q[1];
rz(-1.7439758) q[1];
sx q[1];
rz(1.2417485) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0051796) q[0];
sx q[0];
rz(-1.2997928) q[0];
sx q[0];
rz(-0.47309978) q[0];
rz(-0.21236944) q[2];
sx q[2];
rz(-2.1687733) q[2];
sx q[2];
rz(-0.47456196) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1832402) q[1];
sx q[1];
rz(-1.1179525) q[1];
sx q[1];
rz(-1.1151821) q[1];
x q[2];
rz(-1.0343814) q[3];
sx q[3];
rz(-0.65183833) q[3];
sx q[3];
rz(2.502273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9813098) q[2];
sx q[2];
rz(-0.95429388) q[2];
sx q[2];
rz(2.1295638) q[2];
rz(0.99803287) q[3];
sx q[3];
rz(-0.73603743) q[3];
sx q[3];
rz(1.6691104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0618458) q[0];
sx q[0];
rz(-0.80119067) q[0];
sx q[0];
rz(2.8629942) q[0];
rz(-0.31696907) q[1];
sx q[1];
rz(-1.3895037) q[1];
sx q[1];
rz(0.46151361) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6742548) q[0];
sx q[0];
rz(-1.5227571) q[0];
sx q[0];
rz(-1.4018628) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1912109) q[2];
sx q[2];
rz(-1.912552) q[2];
sx q[2];
rz(2.1555962) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.52109226) q[1];
sx q[1];
rz(-0.25687718) q[1];
sx q[1];
rz(1.9476554) q[1];
x q[2];
rz(-3.0449681) q[3];
sx q[3];
rz(-1.0276252) q[3];
sx q[3];
rz(-0.66978776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3658112) q[2];
sx q[2];
rz(-1.3373988) q[2];
sx q[2];
rz(-2.865045) q[2];
rz(-0.60025275) q[3];
sx q[3];
rz(-0.7163896) q[3];
sx q[3];
rz(-2.6071809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14715956) q[0];
sx q[0];
rz(-0.58013791) q[0];
sx q[0];
rz(-2.4603727) q[0];
rz(-0.45411202) q[1];
sx q[1];
rz(-1.9845767) q[1];
sx q[1];
rz(2.5195794) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2300782) q[0];
sx q[0];
rz(-2.921836) q[0];
sx q[0];
rz(0.4498555) q[0];
rz(0.63779442) q[2];
sx q[2];
rz(-2.0908815) q[2];
sx q[2];
rz(-0.74842036) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4281763) q[1];
sx q[1];
rz(-0.82943664) q[1];
sx q[1];
rz(1.0053289) q[1];
rz(2.8735749) q[3];
sx q[3];
rz(-1.2334196) q[3];
sx q[3];
rz(-1.5806787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6260208) q[2];
sx q[2];
rz(-2.5280648) q[2];
sx q[2];
rz(0.75508368) q[2];
rz(-3.0355022) q[3];
sx q[3];
rz(-1.875149) q[3];
sx q[3];
rz(-2.6809926) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.047091529) q[0];
sx q[0];
rz(-2.5300808) q[0];
sx q[0];
rz(-3.0754572) q[0];
rz(-0.051008929) q[1];
sx q[1];
rz(-2.3936733) q[1];
sx q[1];
rz(1.2640818) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0361657) q[0];
sx q[0];
rz(-1.2107061) q[0];
sx q[0];
rz(2.197916) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35392742) q[2];
sx q[2];
rz(-2.3945237) q[2];
sx q[2];
rz(-0.22176556) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7552774) q[1];
sx q[1];
rz(-1.9615972) q[1];
sx q[1];
rz(0.74797191) q[1];
rz(-pi) q[2];
rz(-2.7080688) q[3];
sx q[3];
rz(-1.2085087) q[3];
sx q[3];
rz(-1.9725245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.85822004) q[2];
sx q[2];
rz(-2.0141352) q[2];
sx q[2];
rz(-0.14632012) q[2];
rz(-0.60394168) q[3];
sx q[3];
rz(-2.5950409) q[3];
sx q[3];
rz(1.7291791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.716575) q[0];
sx q[0];
rz(-1.9442433) q[0];
sx q[0];
rz(-3.0498411) q[0];
rz(3.0692406) q[1];
sx q[1];
rz(-1.3183343) q[1];
sx q[1];
rz(-2.4718557) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7726752) q[0];
sx q[0];
rz(-1.8631422) q[0];
sx q[0];
rz(1.7122853) q[0];
rz(-pi) q[1];
rz(0.40939949) q[2];
sx q[2];
rz(-1.4509307) q[2];
sx q[2];
rz(1.9425962) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.79553855) q[1];
sx q[1];
rz(-0.33949741) q[1];
sx q[1];
rz(0.44793753) q[1];
rz(-pi) q[2];
rz(-2.0132695) q[3];
sx q[3];
rz(-1.4391581) q[3];
sx q[3];
rz(1.7385755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3503795) q[2];
sx q[2];
rz(-2.4243441) q[2];
sx q[2];
rz(-2.4199602) q[2];
rz(-0.77389884) q[3];
sx q[3];
rz(-2.1589203) q[3];
sx q[3];
rz(0.53988808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.72188193) q[0];
sx q[0];
rz(-1.9419436) q[0];
sx q[0];
rz(2.6019959) q[0];
rz(-0.78060141) q[1];
sx q[1];
rz(-1.2465979) q[1];
sx q[1];
rz(0.4221198) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84655428) q[0];
sx q[0];
rz(-1.0031514) q[0];
sx q[0];
rz(-0.66935434) q[0];
x q[1];
rz(-0.49218858) q[2];
sx q[2];
rz(-1.6235311) q[2];
sx q[2];
rz(1.2272138) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.077476689) q[1];
sx q[1];
rz(-1.5327274) q[1];
sx q[1];
rz(-2.0504954) q[1];
x q[2];
rz(2.5529421) q[3];
sx q[3];
rz(-1.6084617) q[3];
sx q[3];
rz(1.879694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0427986) q[2];
sx q[2];
rz(-2.1197987) q[2];
sx q[2];
rz(0.54879028) q[2];
rz(-2.9889066) q[3];
sx q[3];
rz(-2.1137674) q[3];
sx q[3];
rz(0.71167439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3720836) q[0];
sx q[0];
rz(-0.39396572) q[0];
sx q[0];
rz(2.5373996) q[0];
rz(-1.4671885) q[1];
sx q[1];
rz(-1.3507651) q[1];
sx q[1];
rz(-1.1473354) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0096821) q[0];
sx q[0];
rz(-1.1569126) q[0];
sx q[0];
rz(-0.68776806) q[0];
rz(-pi) q[1];
rz(-1.8854972) q[2];
sx q[2];
rz(-2.2813548) q[2];
sx q[2];
rz(1.9942371) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5886138) q[1];
sx q[1];
rz(-1.0880033) q[1];
sx q[1];
rz(1.5881513) q[1];
x q[2];
rz(-0.34965054) q[3];
sx q[3];
rz(-1.9486625) q[3];
sx q[3];
rz(2.1443465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4597822) q[2];
sx q[2];
rz(-1.1638389) q[2];
sx q[2];
rz(-2.5414844) q[2];
rz(-2.4188304) q[3];
sx q[3];
rz(-0.99812752) q[3];
sx q[3];
rz(0.37989894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4031274) q[0];
sx q[0];
rz(-1.5140139) q[0];
sx q[0];
rz(-1.4366666) q[0];
rz(1.5557095) q[1];
sx q[1];
rz(-1.7441505) q[1];
sx q[1];
rz(2.2081262) q[1];
rz(0.49028291) q[2];
sx q[2];
rz(-2.1785707) q[2];
sx q[2];
rz(1.3211801) q[2];
rz(0.71674552) q[3];
sx q[3];
rz(-2.8195753) q[3];
sx q[3];
rz(-1.4618518) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
