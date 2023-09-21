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
rz(2.7058869) q[0];
sx q[0];
rz(11.640179) q[0];
rz(1.9594833) q[1];
sx q[1];
rz(-0.73298454) q[1];
sx q[1];
rz(-2.7690673) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46967888) q[0];
sx q[0];
rz(-2.1436474) q[0];
sx q[0];
rz(-1.8125305) q[0];
rz(-pi) q[1];
rz(-2.9973642) q[2];
sx q[2];
rz(-0.63604522) q[2];
sx q[2];
rz(-0.47338212) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.86245) q[1];
sx q[1];
rz(-1.3819873) q[1];
sx q[1];
rz(0.11629176) q[1];
x q[2];
rz(2.0472237) q[3];
sx q[3];
rz(-2.8989887) q[3];
sx q[3];
rz(-1.1012332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42840502) q[2];
sx q[2];
rz(-1.5593854) q[2];
sx q[2];
rz(0.92450809) q[2];
rz(1.6690147) q[3];
sx q[3];
rz(-1.2481097) q[3];
sx q[3];
rz(1.4991466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18838841) q[0];
sx q[0];
rz(-3.0792455) q[0];
sx q[0];
rz(1.6054608) q[0];
rz(0.19451441) q[1];
sx q[1];
rz(-1.8201927) q[1];
sx q[1];
rz(-0.054873437) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7354436) q[0];
sx q[0];
rz(-0.14649728) q[0];
sx q[0];
rz(-2.2551401) q[0];
rz(0.61972159) q[2];
sx q[2];
rz(-1.1148858) q[2];
sx q[2];
rz(-2.1330619) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4108737) q[1];
sx q[1];
rz(-1.5667856) q[1];
sx q[1];
rz(-1.2282759) q[1];
x q[2];
rz(-1.3589456) q[3];
sx q[3];
rz(-0.90862521) q[3];
sx q[3];
rz(0.3609095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.43869552) q[2];
sx q[2];
rz(-0.56151152) q[2];
sx q[2];
rz(1.4820209) q[2];
rz(-0.99059087) q[3];
sx q[3];
rz(-1.7329268) q[3];
sx q[3];
rz(1.7747442) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3593339) q[0];
sx q[0];
rz(-0.5193091) q[0];
sx q[0];
rz(-2.7666132) q[0];
rz(-2.8938876) q[1];
sx q[1];
rz(-1.9539555) q[1];
sx q[1];
rz(1.8050271) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35555392) q[0];
sx q[0];
rz(-2.9125014) q[0];
sx q[0];
rz(1.8285719) q[0];
rz(-0.61632421) q[2];
sx q[2];
rz(-0.68683544) q[2];
sx q[2];
rz(0.64885215) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3414351) q[1];
sx q[1];
rz(-0.86038024) q[1];
sx q[1];
rz(2.7374949) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7671719) q[3];
sx q[3];
rz(-1.2734405) q[3];
sx q[3];
rz(1.5945827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1027362) q[2];
sx q[2];
rz(-0.83013022) q[2];
sx q[2];
rz(1.0245727) q[2];
rz(-1.2767977) q[3];
sx q[3];
rz(-1.5921311) q[3];
sx q[3];
rz(-1.5156486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17722002) q[0];
sx q[0];
rz(-1.465088) q[0];
sx q[0];
rz(0.033551034) q[0];
rz(-1.230348) q[1];
sx q[1];
rz(-0.76554811) q[1];
sx q[1];
rz(1.4039325) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1407773) q[0];
sx q[0];
rz(-0.40110943) q[0];
sx q[0];
rz(2.6252069) q[0];
rz(-pi) q[1];
rz(0.83909859) q[2];
sx q[2];
rz(-1.027178) q[2];
sx q[2];
rz(2.3630138) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0701323) q[1];
sx q[1];
rz(-1.6825819) q[1];
sx q[1];
rz(1.2791355) q[1];
rz(-pi) q[2];
rz(-2.2622044) q[3];
sx q[3];
rz(-1.5048358) q[3];
sx q[3];
rz(-1.1018167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6491062) q[2];
sx q[2];
rz(-0.84473574) q[2];
sx q[2];
rz(2.9283294) q[2];
rz(3.1212741) q[3];
sx q[3];
rz(-1.3213986) q[3];
sx q[3];
rz(2.7385353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7810818) q[0];
sx q[0];
rz(-1.1163982) q[0];
sx q[0];
rz(-2.1257341) q[0];
rz(1.6332743) q[1];
sx q[1];
rz(-0.95817482) q[1];
sx q[1];
rz(3.1075081) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.236892) q[0];
sx q[0];
rz(-1.1780945) q[0];
sx q[0];
rz(-1.7444201) q[0];
rz(-pi) q[1];
rz(-1.3047406) q[2];
sx q[2];
rz(-1.012946) q[2];
sx q[2];
rz(-1.3706052) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2195697) q[1];
sx q[1];
rz(-1.2864188) q[1];
sx q[1];
rz(-3.1344096) q[1];
rz(-pi) q[2];
rz(2.6328153) q[3];
sx q[3];
rz(-2.1366589) q[3];
sx q[3];
rz(-0.98154991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.74026996) q[2];
sx q[2];
rz(-0.57381845) q[2];
sx q[2];
rz(1.9704698) q[2];
rz(1.3327538) q[3];
sx q[3];
rz(-1.7141902) q[3];
sx q[3];
rz(0.12935054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68833441) q[0];
sx q[0];
rz(-1.1880705) q[0];
sx q[0];
rz(-2.7057498) q[0];
rz(2.0360937) q[1];
sx q[1];
rz(-1.8445797) q[1];
sx q[1];
rz(1.9357392) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3656475) q[0];
sx q[0];
rz(-1.3344814) q[0];
sx q[0];
rz(-2.3274371) q[0];
rz(-pi) q[1];
rz(-1.4756104) q[2];
sx q[2];
rz(-1.5867481) q[2];
sx q[2];
rz(-2.8841281) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1356493) q[1];
sx q[1];
rz(-1.5444396) q[1];
sx q[1];
rz(-0.51075682) q[1];
x q[2];
rz(-1.1711575) q[3];
sx q[3];
rz(-1.6013147) q[3];
sx q[3];
rz(-0.37351028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.897052) q[2];
sx q[2];
rz(-2.5138833) q[2];
sx q[2];
rz(3.0899866) q[2];
rz(0.40766019) q[3];
sx q[3];
rz(-1.9577273) q[3];
sx q[3];
rz(-2.6962962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(2.5200941) q[0];
sx q[0];
rz(-2.5243653) q[0];
sx q[0];
rz(0.67333418) q[0];
rz(0.78397059) q[1];
sx q[1];
rz(-1.4173123) q[1];
sx q[1];
rz(2.6046682) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4250701) q[0];
sx q[0];
rz(-1.5050097) q[0];
sx q[0];
rz(-1.487088) q[0];
rz(-pi) q[1];
rz(-1.3145447) q[2];
sx q[2];
rz(-1.7125687) q[2];
sx q[2];
rz(1.2830551) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.463045) q[1];
sx q[1];
rz(-1.4765413) q[1];
sx q[1];
rz(0.80295347) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7241931) q[3];
sx q[3];
rz(-2.1563357) q[3];
sx q[3];
rz(3.1068387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.15988222) q[2];
sx q[2];
rz(-1.2572224) q[2];
sx q[2];
rz(-0.19101492) q[2];
rz(0.28132176) q[3];
sx q[3];
rz(-1.1912991) q[3];
sx q[3];
rz(-0.34240001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0328338) q[0];
sx q[0];
rz(-2.7712951) q[0];
sx q[0];
rz(-0.38761815) q[0];
rz(-0.11501137) q[1];
sx q[1];
rz(-1.4287881) q[1];
sx q[1];
rz(0.33755916) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4306256) q[0];
sx q[0];
rz(-1.0904878) q[0];
sx q[0];
rz(1.7947012) q[0];
x q[1];
rz(3.0701596) q[2];
sx q[2];
rz(-0.93284235) q[2];
sx q[2];
rz(0.86779867) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3927314) q[1];
sx q[1];
rz(-1.0034605) q[1];
sx q[1];
rz(-0.15565236) q[1];
rz(3.0739325) q[3];
sx q[3];
rz(-2.3849871) q[3];
sx q[3];
rz(0.17186952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.18743029) q[2];
sx q[2];
rz(-2.607441) q[2];
sx q[2];
rz(-1.2672651) q[2];
rz(-0.77504843) q[3];
sx q[3];
rz(-1.5379484) q[3];
sx q[3];
rz(0.78139853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9545814) q[0];
sx q[0];
rz(-1.0382074) q[0];
sx q[0];
rz(0.22098456) q[0];
rz(2.2194608) q[1];
sx q[1];
rz(-1.2698959) q[1];
sx q[1];
rz(-2.1386713) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5054277) q[0];
sx q[0];
rz(-1.822585) q[0];
sx q[0];
rz(-1.7136784) q[0];
rz(2.2641364) q[2];
sx q[2];
rz(-2.4306731) q[2];
sx q[2];
rz(-3.1117709) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3468614) q[1];
sx q[1];
rz(-2.0309629) q[1];
sx q[1];
rz(-0.55333432) q[1];
x q[2];
rz(1.2420197) q[3];
sx q[3];
rz(-2.4283096) q[3];
sx q[3];
rz(2.3438791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6691436) q[2];
sx q[2];
rz(-2.3995235) q[2];
sx q[2];
rz(-2.9830902) q[2];
rz(-2.6878099) q[3];
sx q[3];
rz(-2.3588389) q[3];
sx q[3];
rz(-2.2973072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6145988) q[0];
sx q[0];
rz(-2.232593) q[0];
sx q[0];
rz(2.9504543) q[0];
rz(2.846431) q[1];
sx q[1];
rz(-2.252153) q[1];
sx q[1];
rz(-2.2492762) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.17642) q[0];
sx q[0];
rz(-1.3942379) q[0];
sx q[0];
rz(-2.3202592) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.093776137) q[2];
sx q[2];
rz(-2.2095592) q[2];
sx q[2];
rz(-0.057657777) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.72496966) q[1];
sx q[1];
rz(-1.5048711) q[1];
sx q[1];
rz(-0.23857393) q[1];
rz(-0.51384135) q[3];
sx q[3];
rz(-2.6548879) q[3];
sx q[3];
rz(-2.0972308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7717379) q[2];
sx q[2];
rz(-0.83738804) q[2];
sx q[2];
rz(0.85956335) q[2];
rz(-1.9178948) q[3];
sx q[3];
rz(-0.92646354) q[3];
sx q[3];
rz(-0.74444509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0201465) q[0];
sx q[0];
rz(-2.2059724) q[0];
sx q[0];
rz(-1.7511517) q[0];
rz(1.6620811) q[1];
sx q[1];
rz(-2.6585487) q[1];
sx q[1];
rz(-1.9225635) q[1];
rz(-1.4798726) q[2];
sx q[2];
rz(-0.29279136) q[2];
sx q[2];
rz(1.472483) q[2];
rz(-0.039924351) q[3];
sx q[3];
rz(-1.6394686) q[3];
sx q[3];
rz(-2.6185262) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];