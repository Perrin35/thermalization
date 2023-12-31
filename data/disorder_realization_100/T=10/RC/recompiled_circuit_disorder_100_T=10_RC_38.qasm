OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93958062) q[0];
sx q[0];
rz(-0.35020819) q[0];
sx q[0];
rz(2.7749618) q[0];
rz(-2.2740013) q[1];
sx q[1];
rz(-0.64414135) q[1];
sx q[1];
rz(1.6860513) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6018574) q[0];
sx q[0];
rz(-0.55668133) q[0];
sx q[0];
rz(-2.2599028) q[0];
x q[1];
rz(0.86906616) q[2];
sx q[2];
rz(-1.367374) q[2];
sx q[2];
rz(-0.8806526) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.213233) q[1];
sx q[1];
rz(-1.8395437) q[1];
sx q[1];
rz(-2.1409722) q[1];
x q[2];
rz(1.5815758) q[3];
sx q[3];
rz(-1.618715) q[3];
sx q[3];
rz(2.3673494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.1047487) q[2];
sx q[2];
rz(-1.8061183) q[2];
sx q[2];
rz(2.74995) q[2];
rz(0.13970217) q[3];
sx q[3];
rz(-2.468686) q[3];
sx q[3];
rz(-0.02123775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4166819) q[0];
sx q[0];
rz(-1.4044489) q[0];
sx q[0];
rz(1.2765983) q[0];
rz(-0.87031594) q[1];
sx q[1];
rz(-1.566193) q[1];
sx q[1];
rz(1.2044027) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8613974) q[0];
sx q[0];
rz(-1.9635927) q[0];
sx q[0];
rz(1.1108206) q[0];
rz(0.86512489) q[2];
sx q[2];
rz(-1.7127617) q[2];
sx q[2];
rz(2.1324468) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0871036) q[1];
sx q[1];
rz(-0.37885715) q[1];
sx q[1];
rz(-2.8701251) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9235839) q[3];
sx q[3];
rz(-0.70397607) q[3];
sx q[3];
rz(-2.9670027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5045972) q[2];
sx q[2];
rz(-0.4578788) q[2];
sx q[2];
rz(1.9223928) q[2];
rz(-0.35456625) q[3];
sx q[3];
rz(-0.48186007) q[3];
sx q[3];
rz(-1.569081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5511659) q[0];
sx q[0];
rz(-1.6141491) q[0];
sx q[0];
rz(-0.61808008) q[0];
rz(3.1255426) q[1];
sx q[1];
rz(-2.2647104) q[1];
sx q[1];
rz(1.191167) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1026099) q[0];
sx q[0];
rz(-1.257574) q[0];
sx q[0];
rz(0.15977504) q[0];
x q[1];
rz(-1.3267924) q[2];
sx q[2];
rz(-1.8463496) q[2];
sx q[2];
rz(-0.22305605) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.46398677) q[1];
sx q[1];
rz(-1.7772487) q[1];
sx q[1];
rz(-0.21519214) q[1];
rz(0.54179811) q[3];
sx q[3];
rz(-1.2216611) q[3];
sx q[3];
rz(-2.9475714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1742192) q[2];
sx q[2];
rz(-0.9625532) q[2];
sx q[2];
rz(-1.013914) q[2];
rz(1.3570471) q[3];
sx q[3];
rz(-0.8299399) q[3];
sx q[3];
rz(-1.0323662) q[3];
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
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33092609) q[0];
sx q[0];
rz(-1.7174915) q[0];
sx q[0];
rz(-0.20430918) q[0];
rz(-1.7640242) q[1];
sx q[1];
rz(-1.4148477) q[1];
sx q[1];
rz(-2.2185982) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6072236) q[0];
sx q[0];
rz(-1.2396221) q[0];
sx q[0];
rz(-1.8934728) q[0];
rz(-pi) q[1];
rz(1.5825908) q[2];
sx q[2];
rz(-1.5268541) q[2];
sx q[2];
rz(-0.61461385) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1624565) q[1];
sx q[1];
rz(-2.7297331) q[1];
sx q[1];
rz(-1.3084175) q[1];
rz(-pi) q[2];
rz(-1.5971848) q[3];
sx q[3];
rz(-0.70453405) q[3];
sx q[3];
rz(-1.1772616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4541645) q[2];
sx q[2];
rz(-1.585008) q[2];
sx q[2];
rz(0.50951177) q[2];
rz(2.4984958) q[3];
sx q[3];
rz(-2.0431079) q[3];
sx q[3];
rz(2.174214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61280695) q[0];
sx q[0];
rz(-1.2061773) q[0];
sx q[0];
rz(-2.2221785) q[0];
rz(0.98584229) q[1];
sx q[1];
rz(-1.7835833) q[1];
sx q[1];
rz(1.7808419) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1795792) q[0];
sx q[0];
rz(-0.29173457) q[0];
sx q[0];
rz(-2.4289262) q[0];
rz(-pi) q[1];
rz(-0.64702101) q[2];
sx q[2];
rz(-2.2685452) q[2];
sx q[2];
rz(0.24017142) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6319316) q[1];
sx q[1];
rz(-1.3637378) q[1];
sx q[1];
rz(-0.46405554) q[1];
x q[2];
rz(1.4578044) q[3];
sx q[3];
rz(-0.60243536) q[3];
sx q[3];
rz(3.0756385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.78648606) q[2];
sx q[2];
rz(-1.7559218) q[2];
sx q[2];
rz(3.0267267) q[2];
rz(-2.6857175) q[3];
sx q[3];
rz(-1.3331648) q[3];
sx q[3];
rz(-1.8615287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8961261) q[0];
sx q[0];
rz(-1.2446612) q[0];
sx q[0];
rz(-1.8967569) q[0];
rz(-0.70760977) q[1];
sx q[1];
rz(-1.5039624) q[1];
sx q[1];
rz(-0.27522603) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.595364) q[0];
sx q[0];
rz(-1.9018136) q[0];
sx q[0];
rz(2.0966895) q[0];
rz(-0.92991021) q[2];
sx q[2];
rz(-1.1615796) q[2];
sx q[2];
rz(-2.2355459) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.59906193) q[1];
sx q[1];
rz(-1.9179357) q[1];
sx q[1];
rz(2.6229726) q[1];
rz(-pi) q[2];
rz(-0.082645881) q[3];
sx q[3];
rz(-0.9517037) q[3];
sx q[3];
rz(-1.9649399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.747921) q[2];
sx q[2];
rz(-1.5966406) q[2];
sx q[2];
rz(2.6829524) q[2];
rz(-2.4300872) q[3];
sx q[3];
rz(-0.68370521) q[3];
sx q[3];
rz(0.59035629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8608619) q[0];
sx q[0];
rz(-1.1870528) q[0];
sx q[0];
rz(1.77805) q[0];
rz(-1.4200312) q[1];
sx q[1];
rz(-2.6224711) q[1];
sx q[1];
rz(0.71969676) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2965215) q[0];
sx q[0];
rz(-2.4356027) q[0];
sx q[0];
rz(-0.4476053) q[0];
x q[1];
rz(-2.2816706) q[2];
sx q[2];
rz(-2.2051174) q[2];
sx q[2];
rz(-2.2459522) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4622725) q[1];
sx q[1];
rz(-1.5264891) q[1];
sx q[1];
rz(0.54507749) q[1];
rz(-pi) q[2];
rz(-2.8119874) q[3];
sx q[3];
rz(-0.70809396) q[3];
sx q[3];
rz(-3.0081188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.65934962) q[2];
sx q[2];
rz(-1.6094001) q[2];
sx q[2];
rz(0.68816319) q[2];
rz(-0.041953772) q[3];
sx q[3];
rz(-1.099702) q[3];
sx q[3];
rz(0.28013128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91350895) q[0];
sx q[0];
rz(-1.1627473) q[0];
sx q[0];
rz(-2.9550609) q[0];
rz(-0.60449156) q[1];
sx q[1];
rz(-1.0124606) q[1];
sx q[1];
rz(-1.6465181) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5160617) q[0];
sx q[0];
rz(-0.4058668) q[0];
sx q[0];
rz(-0.56674515) q[0];
x q[1];
rz(0.69627701) q[2];
sx q[2];
rz(-1.7211203) q[2];
sx q[2];
rz(2.6339649) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9547792) q[1];
sx q[1];
rz(-2.8615132) q[1];
sx q[1];
rz(-0.38444744) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64077326) q[3];
sx q[3];
rz(-0.91859222) q[3];
sx q[3];
rz(-0.81354248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8994393) q[2];
sx q[2];
rz(-2.184325) q[2];
sx q[2];
rz(0.014766679) q[2];
rz(-1.2773369) q[3];
sx q[3];
rz(-1.3093964) q[3];
sx q[3];
rz(1.4130672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9799141) q[0];
sx q[0];
rz(-2.4665687) q[0];
sx q[0];
rz(0.87798464) q[0];
rz(-0.19628482) q[1];
sx q[1];
rz(-1.1801964) q[1];
sx q[1];
rz(1.3605798) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.774705) q[0];
sx q[0];
rz(-0.57292143) q[0];
sx q[0];
rz(1.6402871) q[0];
rz(-pi) q[1];
rz(-2.1979245) q[2];
sx q[2];
rz(-1.6086279) q[2];
sx q[2];
rz(0.19247069) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0226375) q[1];
sx q[1];
rz(-1.9070909) q[1];
sx q[1];
rz(2.789546) q[1];
rz(0.15501539) q[3];
sx q[3];
rz(-1.1559556) q[3];
sx q[3];
rz(2.993194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.62375162) q[2];
sx q[2];
rz(-1.6343642) q[2];
sx q[2];
rz(0.8927792) q[2];
rz(0.039285224) q[3];
sx q[3];
rz(-1.4792484) q[3];
sx q[3];
rz(-2.3915496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84353012) q[0];
sx q[0];
rz(-3.1382882) q[0];
sx q[0];
rz(3.0950586) q[0];
rz(2.2325366) q[1];
sx q[1];
rz(-0.87288705) q[1];
sx q[1];
rz(0.7199026) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4601446) q[0];
sx q[0];
rz(-1.3272734) q[0];
sx q[0];
rz(0.081865099) q[0];
x q[1];
rz(-1.5653531) q[2];
sx q[2];
rz(-1.8796225) q[2];
sx q[2];
rz(0.18143166) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.1137005) q[1];
sx q[1];
rz(-2.8942331) q[1];
sx q[1];
rz(-1.6340096) q[1];
rz(0.73086892) q[3];
sx q[3];
rz(-2.1248397) q[3];
sx q[3];
rz(-2.4027367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4204734) q[2];
sx q[2];
rz(-1.7420008) q[2];
sx q[2];
rz(-3.1151248) q[2];
rz(1.8376393) q[3];
sx q[3];
rz(-2.3243258) q[3];
sx q[3];
rz(1.8855689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7883041) q[0];
sx q[0];
rz(-1.7699387) q[0];
sx q[0];
rz(1.8040245) q[0];
rz(2.9251255) q[1];
sx q[1];
rz(-1.6995866) q[1];
sx q[1];
rz(1.235984) q[1];
rz(-1.886006) q[2];
sx q[2];
rz(-0.8885347) q[2];
sx q[2];
rz(2.3226429) q[2];
rz(-2.1574216) q[3];
sx q[3];
rz(-1.5979206) q[3];
sx q[3];
rz(0.4978705) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
