OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.714158) q[0];
sx q[0];
rz(-2.7058869) q[0];
sx q[0];
rz(-0.92619196) q[0];
rz(-1.1821094) q[1];
sx q[1];
rz(-2.4086081) q[1];
sx q[1];
rz(2.7690673) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0988136) q[0];
sx q[0];
rz(-2.5251237) q[0];
sx q[0];
rz(2.7862076) q[0];
rz(-pi) q[1];
rz(2.9973642) q[2];
sx q[2];
rz(-0.63604522) q[2];
sx q[2];
rz(0.47338212) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2791427) q[1];
sx q[1];
rz(-1.3819873) q[1];
sx q[1];
rz(0.11629176) q[1];
x q[2];
rz(-3.0285809) q[3];
sx q[3];
rz(-1.7859308) q[3];
sx q[3];
rz(-2.5290031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42840502) q[2];
sx q[2];
rz(-1.5593854) q[2];
sx q[2];
rz(-2.2170846) q[2];
rz(1.472578) q[3];
sx q[3];
rz(-1.893483) q[3];
sx q[3];
rz(1.4991466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18838841) q[0];
sx q[0];
rz(-0.062347118) q[0];
sx q[0];
rz(-1.5361319) q[0];
rz(2.9470782) q[1];
sx q[1];
rz(-1.8201927) q[1];
sx q[1];
rz(0.054873437) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48557278) q[0];
sx q[0];
rz(-1.6632073) q[0];
sx q[0];
rz(1.4569605) q[0];
x q[1];
rz(-2.4404281) q[2];
sx q[2];
rz(-0.75116457) q[2];
sx q[2];
rz(3.1322111) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.730719) q[1];
sx q[1];
rz(-1.5667856) q[1];
sx q[1];
rz(1.9133168) q[1];
x q[2];
rz(-2.4684286) q[3];
sx q[3];
rz(-1.737397) q[3];
sx q[3];
rz(1.8002321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7028971) q[2];
sx q[2];
rz(-2.5800811) q[2];
sx q[2];
rz(-1.4820209) q[2];
rz(2.1510018) q[3];
sx q[3];
rz(-1.4086658) q[3];
sx q[3];
rz(-1.7747442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7822587) q[0];
sx q[0];
rz(-0.5193091) q[0];
sx q[0];
rz(0.3749795) q[0];
rz(-2.8938876) q[1];
sx q[1];
rz(-1.9539555) q[1];
sx q[1];
rz(-1.3365655) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5216615) q[0];
sx q[0];
rz(-1.3494028) q[0];
sx q[0];
rz(-3.0822166) q[0];
x q[1];
rz(-0.61632421) q[2];
sx q[2];
rz(-2.4547572) q[2];
sx q[2];
rz(2.4927405) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6429236) q[1];
sx q[1];
rz(-1.8734697) q[1];
sx q[1];
rz(-2.32294) q[1];
x q[2];
rz(-0.7671719) q[3];
sx q[3];
rz(-1.8681521) q[3];
sx q[3];
rz(-1.5470099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1027362) q[2];
sx q[2];
rz(-2.3114624) q[2];
sx q[2];
rz(-2.1170199) q[2];
rz(-1.2767977) q[3];
sx q[3];
rz(-1.5921311) q[3];
sx q[3];
rz(1.6259441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9643726) q[0];
sx q[0];
rz(-1.465088) q[0];
sx q[0];
rz(-3.1080416) q[0];
rz(-1.230348) q[1];
sx q[1];
rz(-0.76554811) q[1];
sx q[1];
rz(1.4039325) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1407773) q[0];
sx q[0];
rz(-2.7404832) q[0];
sx q[0];
rz(0.51638575) q[0];
x q[1];
rz(-2.4594175) q[2];
sx q[2];
rz(-0.96207679) q[2];
sx q[2];
rz(-1.9145554) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.85642203) q[1];
sx q[1];
rz(-0.31177786) q[1];
sx q[1];
rz(1.1986033) q[1];
x q[2];
rz(-2.2622044) q[3];
sx q[3];
rz(-1.5048358) q[3];
sx q[3];
rz(2.0397759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6491062) q[2];
sx q[2];
rz(-0.84473574) q[2];
sx q[2];
rz(0.21326324) q[2];
rz(3.1212741) q[3];
sx q[3];
rz(-1.3213986) q[3];
sx q[3];
rz(-0.4030574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3605109) q[0];
sx q[0];
rz(-1.1163982) q[0];
sx q[0];
rz(1.0158585) q[0];
rz(-1.6332743) q[1];
sx q[1];
rz(-2.1834178) q[1];
sx q[1];
rz(-0.034084592) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7406697) q[0];
sx q[0];
rz(-1.4105083) q[0];
sx q[0];
rz(0.39808654) q[0];
x q[1];
rz(-2.5675315) q[2];
sx q[2];
rz(-1.3458369) q[2];
sx q[2];
rz(-0.056919295) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7923813) q[1];
sx q[1];
rz(-1.5639018) q[1];
sx q[1];
rz(-1.8551808) q[1];
x q[2];
rz(-2.2250416) q[3];
sx q[3];
rz(-2.399728) q[3];
sx q[3];
rz(-1.3548917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4013227) q[2];
sx q[2];
rz(-2.5677742) q[2];
sx q[2];
rz(-1.1711228) q[2];
rz(1.8088388) q[3];
sx q[3];
rz(-1.7141902) q[3];
sx q[3];
rz(3.0122421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4532582) q[0];
sx q[0];
rz(-1.1880705) q[0];
sx q[0];
rz(-2.7057498) q[0];
rz(1.1054989) q[1];
sx q[1];
rz(-1.8445797) q[1];
sx q[1];
rz(-1.9357392) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0379369) q[0];
sx q[0];
rz(-0.78563848) q[0];
sx q[0];
rz(1.908179) q[0];
x q[1];
rz(3.1255683) q[2];
sx q[2];
rz(-1.4756225) q[2];
sx q[2];
rz(1.3118088) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7534415) q[1];
sx q[1];
rz(-2.6302164) q[1];
sx q[1];
rz(-0.053877342) q[1];
x q[2];
rz(-1.6490963) q[3];
sx q[3];
rz(-0.40073985) q[3];
sx q[3];
rz(-2.0164255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.897052) q[2];
sx q[2];
rz(-0.62770939) q[2];
sx q[2];
rz(-0.051606027) q[2];
rz(0.40766019) q[3];
sx q[3];
rz(-1.1838653) q[3];
sx q[3];
rz(2.6962962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5200941) q[0];
sx q[0];
rz(-2.5243653) q[0];
sx q[0];
rz(2.4682585) q[0];
rz(-0.78397059) q[1];
sx q[1];
rz(-1.4173123) q[1];
sx q[1];
rz(-2.6046682) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4250701) q[0];
sx q[0];
rz(-1.5050097) q[0];
sx q[0];
rz(-1.487088) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3145447) q[2];
sx q[2];
rz(-1.7125687) q[2];
sx q[2];
rz(-1.8585376) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6785477) q[1];
sx q[1];
rz(-1.4765413) q[1];
sx q[1];
rz(0.80295347) q[1];
rz(1.4173996) q[3];
sx q[3];
rz(-0.985257) q[3];
sx q[3];
rz(-0.034754001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9817104) q[2];
sx q[2];
rz(-1.8843702) q[2];
sx q[2];
rz(-0.19101492) q[2];
rz(0.28132176) q[3];
sx q[3];
rz(-1.9502935) q[3];
sx q[3];
rz(-2.7991926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0328338) q[0];
sx q[0];
rz(-0.37029752) q[0];
sx q[0];
rz(2.7539745) q[0];
rz(3.0265813) q[1];
sx q[1];
rz(-1.4287881) q[1];
sx q[1];
rz(-2.8040335) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8885376) q[0];
sx q[0];
rz(-2.6153784) q[0];
sx q[0];
rz(-2.7387268) q[0];
x q[1];
rz(-2.2099724) q[2];
sx q[2];
rz(-1.5134303) q[2];
sx q[2];
rz(0.66040874) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3927314) q[1];
sx q[1];
rz(-1.0034605) q[1];
sx q[1];
rz(0.15565236) q[1];
x q[2];
rz(1.5070595) q[3];
sx q[3];
rz(-2.3252441) q[3];
sx q[3];
rz(3.0626429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9541624) q[2];
sx q[2];
rz(-0.53415161) q[2];
sx q[2];
rz(-1.2672651) q[2];
rz(-2.3665442) q[3];
sx q[3];
rz(-1.5379484) q[3];
sx q[3];
rz(2.3601941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18701126) q[0];
sx q[0];
rz(-2.1033852) q[0];
sx q[0];
rz(-2.9206081) q[0];
rz(0.9221319) q[1];
sx q[1];
rz(-1.8716967) q[1];
sx q[1];
rz(-2.1386713) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1598338) q[0];
sx q[0];
rz(-2.8528385) q[0];
sx q[0];
rz(-2.6360378) q[0];
x q[1];
rz(-0.9858176) q[2];
sx q[2];
rz(-1.1406116) q[2];
sx q[2];
rz(-2.1625724) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.79473125) q[1];
sx q[1];
rz(-2.0309629) q[1];
sx q[1];
rz(-0.55333432) q[1];
rz(-pi) q[2];
rz(2.2569611) q[3];
sx q[3];
rz(-1.7836708) q[3];
sx q[3];
rz(-0.5205982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4724491) q[2];
sx q[2];
rz(-0.74206918) q[2];
sx q[2];
rz(0.15850244) q[2];
rz(-0.45378271) q[3];
sx q[3];
rz(-2.3588389) q[3];
sx q[3];
rz(2.2973072) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52699387) q[0];
sx q[0];
rz(-0.90899962) q[0];
sx q[0];
rz(0.19113834) q[0];
rz(-2.846431) q[1];
sx q[1];
rz(-0.8894397) q[1];
sx q[1];
rz(-2.2492762) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3494209) q[0];
sx q[0];
rz(-0.76602174) q[0];
sx q[0];
rz(-1.3146521) q[0];
rz(-pi) q[1];
rz(3.0478165) q[2];
sx q[2];
rz(-2.2095592) q[2];
sx q[2];
rz(3.0839349) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.72496966) q[1];
sx q[1];
rz(-1.6367216) q[1];
sx q[1];
rz(-0.23857393) q[1];
x q[2];
rz(1.3163371) q[3];
sx q[3];
rz(-1.1512827) q[3];
sx q[3];
rz(-2.6655243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7717379) q[2];
sx q[2];
rz(-0.83738804) q[2];
sx q[2];
rz(-0.85956335) q[2];
rz(1.2236979) q[3];
sx q[3];
rz(-0.92646354) q[3];
sx q[3];
rz(-0.74444509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1214462) q[0];
sx q[0];
rz(-0.93562026) q[0];
sx q[0];
rz(1.390441) q[0];
rz(1.4795115) q[1];
sx q[1];
rz(-0.48304396) q[1];
sx q[1];
rz(1.2190291) q[1];
rz(1.4798726) q[2];
sx q[2];
rz(-2.8488013) q[2];
sx q[2];
rz(-1.6691096) q[2];
rz(-2.0966093) q[3];
sx q[3];
rz(-0.079418728) q[3];
sx q[3];
rz(3.1374745) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
