OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4274347) q[0];
sx q[0];
rz(-0.43570575) q[0];
sx q[0];
rz(-2.2154007) q[0];
rz(1.9594833) q[1];
sx q[1];
rz(-0.73298454) q[1];
sx q[1];
rz(0.37252537) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0988136) q[0];
sx q[0];
rz(-0.61646898) q[0];
sx q[0];
rz(-2.7862076) q[0];
rz(1.6765321) q[2];
sx q[2];
rz(-2.1991962) q[2];
sx q[2];
rz(2.489593) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.86245) q[1];
sx q[1];
rz(-1.3819873) q[1];
sx q[1];
rz(-3.0253009) q[1];
rz(-pi) q[2];
x q[2];
rz(1.094369) q[3];
sx q[3];
rz(-0.24260394) q[3];
sx q[3];
rz(2.0403595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.42840502) q[2];
sx q[2];
rz(-1.5593854) q[2];
sx q[2];
rz(-2.2170846) q[2];
rz(-1.472578) q[3];
sx q[3];
rz(-1.893483) q[3];
sx q[3];
rz(1.6424461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18838841) q[0];
sx q[0];
rz(-3.0792455) q[0];
sx q[0];
rz(-1.5361319) q[0];
rz(0.19451441) q[1];
sx q[1];
rz(-1.3214) q[1];
sx q[1];
rz(-3.0867192) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40614906) q[0];
sx q[0];
rz(-2.9950954) q[0];
sx q[0];
rz(-2.2551401) q[0];
rz(-pi) q[1];
rz(2.4404281) q[2];
sx q[2];
rz(-0.75116457) q[2];
sx q[2];
rz(-3.1322111) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4108737) q[1];
sx q[1];
rz(-1.5748071) q[1];
sx q[1];
rz(1.2282759) q[1];
rz(-pi) q[2];
rz(-0.67316405) q[3];
sx q[3];
rz(-1.737397) q[3];
sx q[3];
rz(1.3413606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43869552) q[2];
sx q[2];
rz(-2.5800811) q[2];
sx q[2];
rz(-1.6595718) q[2];
rz(-2.1510018) q[3];
sx q[3];
rz(-1.7329268) q[3];
sx q[3];
rz(-1.7747442) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7822587) q[0];
sx q[0];
rz(-2.6222836) q[0];
sx q[0];
rz(0.3749795) q[0];
rz(-2.8938876) q[1];
sx q[1];
rz(-1.9539555) q[1];
sx q[1];
rz(-1.3365655) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96391812) q[0];
sx q[0];
rz(-1.5128711) q[0];
sx q[0];
rz(1.3490246) q[0];
rz(-pi) q[1];
rz(2.5518718) q[2];
sx q[2];
rz(-1.1955185) q[2];
sx q[2];
rz(-0.42082618) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.49866906) q[1];
sx q[1];
rz(-1.2681229) q[1];
sx q[1];
rz(-0.81865262) q[1];
rz(2.3744208) q[3];
sx q[3];
rz(-1.8681521) q[3];
sx q[3];
rz(1.5945827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0388564) q[2];
sx q[2];
rz(-0.83013022) q[2];
sx q[2];
rz(-1.0245727) q[2];
rz(1.864795) q[3];
sx q[3];
rz(-1.5921311) q[3];
sx q[3];
rz(1.6259441) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9643726) q[0];
sx q[0];
rz(-1.6765046) q[0];
sx q[0];
rz(-0.033551034) q[0];
rz(1.230348) q[1];
sx q[1];
rz(-2.3760445) q[1];
sx q[1];
rz(1.4039325) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0899635) q[0];
sx q[0];
rz(-1.7647867) q[0];
sx q[0];
rz(2.7882663) q[0];
x q[1];
rz(-2.3061182) q[2];
sx q[2];
rz(-0.88047853) q[2];
sx q[2];
rz(-2.8719605) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0701323) q[1];
sx q[1];
rz(-1.4590108) q[1];
sx q[1];
rz(-1.8624572) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2622044) q[3];
sx q[3];
rz(-1.5048358) q[3];
sx q[3];
rz(-1.1018167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4924865) q[2];
sx q[2];
rz(-2.2968569) q[2];
sx q[2];
rz(2.9283294) q[2];
rz(3.1212741) q[3];
sx q[3];
rz(-1.820194) q[3];
sx q[3];
rz(-2.7385353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3605109) q[0];
sx q[0];
rz(-1.1163982) q[0];
sx q[0];
rz(-1.0158585) q[0];
rz(-1.5083183) q[1];
sx q[1];
rz(-0.95817482) q[1];
sx q[1];
rz(3.1075081) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.236892) q[0];
sx q[0];
rz(-1.9634982) q[0];
sx q[0];
rz(1.3971726) q[0];
rz(-1.3047406) q[2];
sx q[2];
rz(-2.1286466) q[2];
sx q[2];
rz(-1.7709874) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8964256) q[1];
sx q[1];
rz(-0.28446576) q[1];
sx q[1];
rz(1.5953654) q[1];
rz(0.94200763) q[3];
sx q[3];
rz(-1.1470456) q[3];
sx q[3];
rz(-2.8429192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4013227) q[2];
sx q[2];
rz(-2.5677742) q[2];
sx q[2];
rz(1.9704698) q[2];
rz(-1.8088388) q[3];
sx q[3];
rz(-1.4274024) q[3];
sx q[3];
rz(3.0122421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68833441) q[0];
sx q[0];
rz(-1.1880705) q[0];
sx q[0];
rz(2.7057498) q[0];
rz(-2.0360937) q[1];
sx q[1];
rz(-1.8445797) q[1];
sx q[1];
rz(1.2058535) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1036557) q[0];
sx q[0];
rz(-2.3559542) q[0];
sx q[0];
rz(1.908179) q[0];
rz(-pi) q[1];
rz(1.4756104) q[2];
sx q[2];
rz(-1.5867481) q[2];
sx q[2];
rz(-0.25746458) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7534415) q[1];
sx q[1];
rz(-2.6302164) q[1];
sx q[1];
rz(0.053877342) q[1];
rz(-pi) q[2];
rz(-1.6490963) q[3];
sx q[3];
rz(-2.7408528) q[3];
sx q[3];
rz(-1.1251671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2445406) q[2];
sx q[2];
rz(-2.5138833) q[2];
sx q[2];
rz(-3.0899866) q[2];
rz(-0.40766019) q[3];
sx q[3];
rz(-1.1838653) q[3];
sx q[3];
rz(-2.6962962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5200941) q[0];
sx q[0];
rz(-2.5243653) q[0];
sx q[0];
rz(2.4682585) q[0];
rz(2.3576221) q[1];
sx q[1];
rz(-1.7242804) q[1];
sx q[1];
rz(-0.53692445) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2818031) q[0];
sx q[0];
rz(-1.6543232) q[0];
sx q[0];
rz(-0.066017166) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9951018) q[2];
sx q[2];
rz(-1.8244201) q[2];
sx q[2];
rz(-0.25073642) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.19837025) q[1];
sx q[1];
rz(-2.3343625) q[1];
sx q[1];
rz(0.13065773) q[1];
x q[2];
rz(-2.9151239) q[3];
sx q[3];
rz(-0.60301757) q[3];
sx q[3];
rz(-2.9035567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9817104) q[2];
sx q[2];
rz(-1.8843702) q[2];
sx q[2];
rz(-2.9505777) q[2];
rz(-2.8602709) q[3];
sx q[3];
rz(-1.9502935) q[3];
sx q[3];
rz(0.34240001) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1087588) q[0];
sx q[0];
rz(-0.37029752) q[0];
sx q[0];
rz(-2.7539745) q[0];
rz(-3.0265813) q[1];
sx q[1];
rz(-1.4287881) q[1];
sx q[1];
rz(2.8040335) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8885376) q[0];
sx q[0];
rz(-0.5262143) q[0];
sx q[0];
rz(2.7387268) q[0];
rz(-0.93162025) q[2];
sx q[2];
rz(-1.5134303) q[2];
sx q[2];
rz(-0.66040874) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2355289) q[1];
sx q[1];
rz(-1.7019094) q[1];
sx q[1];
rz(-2.1436585) q[1];
rz(-pi) q[2];
x q[2];
rz(0.067660178) q[3];
sx q[3];
rz(-2.3849871) q[3];
sx q[3];
rz(2.9697231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9541624) q[2];
sx q[2];
rz(-2.607441) q[2];
sx q[2];
rz(-1.2672651) q[2];
rz(0.77504843) q[3];
sx q[3];
rz(-1.5379484) q[3];
sx q[3];
rz(2.3601941) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9545814) q[0];
sx q[0];
rz(-2.1033852) q[0];
sx q[0];
rz(2.9206081) q[0];
rz(-0.9221319) q[1];
sx q[1];
rz(-1.2698959) q[1];
sx q[1];
rz(-2.1386713) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.040398) q[0];
sx q[0];
rz(-1.4324491) q[0];
sx q[0];
rz(0.25427108) q[0];
x q[1];
rz(-2.6384764) q[2];
sx q[2];
rz(-2.0965577) q[2];
sx q[2];
rz(-0.8612649) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3468614) q[1];
sx q[1];
rz(-1.1106297) q[1];
sx q[1];
rz(0.55333432) q[1];
rz(-pi) q[2];
rz(2.2569611) q[3];
sx q[3];
rz(-1.7836708) q[3];
sx q[3];
rz(2.6209944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6691436) q[2];
sx q[2];
rz(-2.3995235) q[2];
sx q[2];
rz(2.9830902) q[2];
rz(2.6878099) q[3];
sx q[3];
rz(-2.3588389) q[3];
sx q[3];
rz(-0.84428549) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52699387) q[0];
sx q[0];
rz(-2.232593) q[0];
sx q[0];
rz(0.19113834) q[0];
rz(2.846431) q[1];
sx q[1];
rz(-2.252153) q[1];
sx q[1];
rz(-2.2492762) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79217171) q[0];
sx q[0];
rz(-0.76602174) q[0];
sx q[0];
rz(-1.3146521) q[0];
x q[1];
rz(0.92992444) q[2];
sx q[2];
rz(-1.495549) q[2];
sx q[2];
rz(-1.6844695) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.72496966) q[1];
sx q[1];
rz(-1.6367216) q[1];
sx q[1];
rz(-2.9030187) q[1];
rz(-pi) q[2];
rz(-1.8252556) q[3];
sx q[3];
rz(-1.99031) q[3];
sx q[3];
rz(-0.47606836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3698547) q[2];
sx q[2];
rz(-2.3042046) q[2];
sx q[2];
rz(-2.2820293) q[2];
rz(-1.2236979) q[3];
sx q[3];
rz(-0.92646354) q[3];
sx q[3];
rz(0.74444509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1214462) q[0];
sx q[0];
rz(-0.93562026) q[0];
sx q[0];
rz(1.390441) q[0];
rz(-1.4795115) q[1];
sx q[1];
rz(-2.6585487) q[1];
sx q[1];
rz(-1.9225635) q[1];
rz(1.2791469) q[2];
sx q[2];
rz(-1.5970061) q[2];
sx q[2];
rz(2.9562052) q[2];
rz(-3.1016683) q[3];
sx q[3];
rz(-1.502124) q[3];
sx q[3];
rz(0.52306642) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
