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
rz(-2.7690673) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2339708) q[0];
sx q[0];
rz(-1.3682433) q[0];
sx q[0];
rz(-2.5552208) q[0];
rz(-pi) q[1];
rz(-1.4650605) q[2];
sx q[2];
rz(-2.1991962) q[2];
sx q[2];
rz(-0.65199967) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.86245) q[1];
sx q[1];
rz(-1.3819873) q[1];
sx q[1];
rz(0.11629176) q[1];
rz(1.7872693) q[3];
sx q[3];
rz(-1.6811922) q[3];
sx q[3];
rz(0.9339827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7131876) q[2];
sx q[2];
rz(-1.5593854) q[2];
sx q[2];
rz(2.2170846) q[2];
rz(1.6690147) q[3];
sx q[3];
rz(-1.2481097) q[3];
sx q[3];
rz(-1.6424461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(0.18838841) q[0];
sx q[0];
rz(-0.062347118) q[0];
sx q[0];
rz(-1.5361319) q[0];
rz(2.9470782) q[1];
sx q[1];
rz(-1.8201927) q[1];
sx q[1];
rz(-3.0867192) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0458192) q[0];
sx q[0];
rz(-1.4574483) q[0];
sx q[0];
rz(-0.093009526) q[0];
rz(1.0286249) q[2];
sx q[2];
rz(-1.0222058) q[2];
sx q[2];
rz(-0.86663914) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3127628) q[1];
sx q[1];
rz(-0.34254303) q[1];
sx q[1];
rz(-1.5827375) q[1];
rz(-0.67316405) q[3];
sx q[3];
rz(-1.737397) q[3];
sx q[3];
rz(1.3413606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7028971) q[2];
sx q[2];
rz(-0.56151152) q[2];
sx q[2];
rz(-1.6595718) q[2];
rz(-2.1510018) q[3];
sx q[3];
rz(-1.4086658) q[3];
sx q[3];
rz(-1.3668485) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3593339) q[0];
sx q[0];
rz(-0.5193091) q[0];
sx q[0];
rz(2.7666132) q[0];
rz(2.8938876) q[1];
sx q[1];
rz(-1.9539555) q[1];
sx q[1];
rz(1.3365655) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96391812) q[0];
sx q[0];
rz(-1.6287215) q[0];
sx q[0];
rz(-1.3490246) q[0];
x q[1];
rz(0.61632421) q[2];
sx q[2];
rz(-2.4547572) q[2];
sx q[2];
rz(-2.4927405) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.8001576) q[1];
sx q[1];
rz(-0.86038024) q[1];
sx q[1];
rz(2.7374949) q[1];
rz(-pi) q[2];
rz(-1.1683447) q[3];
sx q[3];
rz(-0.84512049) q[3];
sx q[3];
rz(0.25154164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0388564) q[2];
sx q[2];
rz(-0.83013022) q[2];
sx q[2];
rz(-2.1170199) q[2];
rz(-1.864795) q[3];
sx q[3];
rz(-1.5921311) q[3];
sx q[3];
rz(1.5156486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17722002) q[0];
sx q[0];
rz(-1.465088) q[0];
sx q[0];
rz(-0.033551034) q[0];
rz(-1.230348) q[1];
sx q[1];
rz(-2.3760445) q[1];
sx q[1];
rz(1.7376602) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1407773) q[0];
sx q[0];
rz(-2.7404832) q[0];
sx q[0];
rz(0.51638575) q[0];
rz(-2.3061182) q[2];
sx q[2];
rz(-2.2611141) q[2];
sx q[2];
rz(-0.26963216) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.4671859) q[1];
sx q[1];
rz(-1.8605839) q[1];
sx q[1];
rz(0.11667103) q[1];
x q[2];
rz(3.0560533) q[3];
sx q[3];
rz(-0.88118689) q[3];
sx q[3];
rz(-0.52348189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4924865) q[2];
sx q[2];
rz(-0.84473574) q[2];
sx q[2];
rz(0.21326324) q[2];
rz(-3.1212741) q[3];
sx q[3];
rz(-1.3213986) q[3];
sx q[3];
rz(0.4030574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3605109) q[0];
sx q[0];
rz(-1.1163982) q[0];
sx q[0];
rz(-1.0158585) q[0];
rz(1.6332743) q[1];
sx q[1];
rz(-0.95817482) q[1];
sx q[1];
rz(3.1075081) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7406697) q[0];
sx q[0];
rz(-1.4105083) q[0];
sx q[0];
rz(-0.39808654) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7427865) q[2];
sx q[2];
rz(-2.5296693) q[2];
sx q[2];
rz(-1.8460225) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9220229) q[1];
sx q[1];
rz(-1.2864188) q[1];
sx q[1];
rz(3.1344096) q[1];
rz(-0.91655101) q[3];
sx q[3];
rz(-0.74186462) q[3];
sx q[3];
rz(1.786701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.74026996) q[2];
sx q[2];
rz(-0.57381845) q[2];
sx q[2];
rz(1.1711228) q[2];
rz(1.8088388) q[3];
sx q[3];
rz(-1.7141902) q[3];
sx q[3];
rz(-0.12935054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4532582) q[0];
sx q[0];
rz(-1.9535221) q[0];
sx q[0];
rz(-0.43584287) q[0];
rz(2.0360937) q[1];
sx q[1];
rz(-1.8445797) q[1];
sx q[1];
rz(-1.2058535) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5642731) q[0];
sx q[0];
rz(-2.3015129) q[0];
sx q[0];
rz(0.3198091) q[0];
rz(-pi) q[1];
rz(-1.7370991) q[2];
sx q[2];
rz(-0.096509343) q[2];
sx q[2];
rz(1.6627179) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.3881512) q[1];
sx q[1];
rz(-2.6302164) q[1];
sx q[1];
rz(-3.0877153) q[1];
rz(-1.9704351) q[3];
sx q[3];
rz(-1.6013147) q[3];
sx q[3];
rz(0.37351028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2445406) q[2];
sx q[2];
rz(-2.5138833) q[2];
sx q[2];
rz(-0.051606027) q[2];
rz(-0.40766019) q[3];
sx q[3];
rz(-1.1838653) q[3];
sx q[3];
rz(-2.6962962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62149858) q[0];
sx q[0];
rz(-0.61722732) q[0];
sx q[0];
rz(-2.4682585) q[0];
rz(0.78397059) q[1];
sx q[1];
rz(-1.7242804) q[1];
sx q[1];
rz(-2.6046682) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4250701) q[0];
sx q[0];
rz(-1.636583) q[0];
sx q[0];
rz(-1.487088) q[0];
x q[1];
rz(-1.3145447) q[2];
sx q[2];
rz(-1.429024) q[2];
sx q[2];
rz(1.8585376) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.19837025) q[1];
sx q[1];
rz(-0.80723019) q[1];
sx q[1];
rz(3.0109349) q[1];
rz(-pi) q[2];
rz(1.7241931) q[3];
sx q[3];
rz(-0.985257) q[3];
sx q[3];
rz(-3.1068387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.15988222) q[2];
sx q[2];
rz(-1.8843702) q[2];
sx q[2];
rz(2.9505777) q[2];
rz(-0.28132176) q[3];
sx q[3];
rz(-1.9502935) q[3];
sx q[3];
rz(2.7991926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1087588) q[0];
sx q[0];
rz(-2.7712951) q[0];
sx q[0];
rz(2.7539745) q[0];
rz(-0.11501137) q[1];
sx q[1];
rz(-1.4287881) q[1];
sx q[1];
rz(-2.8040335) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4306256) q[0];
sx q[0];
rz(-2.0511048) q[0];
sx q[0];
rz(1.3468915) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6667716) q[2];
sx q[2];
rz(-0.64138597) q[2];
sx q[2];
rz(-0.98737398) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.464752) q[1];
sx q[1];
rz(-0.58603474) q[1];
sx q[1];
rz(1.3321484) q[1];
rz(-1.6345331) q[3];
sx q[3];
rz(-2.3252441) q[3];
sx q[3];
rz(3.0626429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.18743029) q[2];
sx q[2];
rz(-0.53415161) q[2];
sx q[2];
rz(-1.8743275) q[2];
rz(0.77504843) q[3];
sx q[3];
rz(-1.6036443) q[3];
sx q[3];
rz(-2.3601941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.18701126) q[0];
sx q[0];
rz(-1.0382074) q[0];
sx q[0];
rz(2.9206081) q[0];
rz(-2.2194608) q[1];
sx q[1];
rz(-1.2698959) q[1];
sx q[1];
rz(2.1386713) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1598338) q[0];
sx q[0];
rz(-2.8528385) q[0];
sx q[0];
rz(-2.6360378) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1557751) q[2];
sx q[2];
rz(-1.1406116) q[2];
sx q[2];
rz(-0.9790203) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6332615) q[1];
sx q[1];
rz(-1.08053) q[1];
sx q[1];
rz(-2.0983178) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27243154) q[3];
sx q[3];
rz(-2.2386132) q[3];
sx q[3];
rz(1.2215134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4724491) q[2];
sx q[2];
rz(-2.3995235) q[2];
sx q[2];
rz(-0.15850244) q[2];
rz(2.6878099) q[3];
sx q[3];
rz(-2.3588389) q[3];
sx q[3];
rz(-0.84428549) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6145988) q[0];
sx q[0];
rz(-0.90899962) q[0];
sx q[0];
rz(0.19113834) q[0];
rz(-0.29516164) q[1];
sx q[1];
rz(-0.8894397) q[1];
sx q[1];
rz(2.2492762) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6979881) q[0];
sx q[0];
rz(-0.83570489) q[0];
sx q[0];
rz(-2.9025335) q[0];
rz(1.4453663) q[2];
sx q[2];
rz(-2.4969366) q[2];
sx q[2];
rz(0.21412011) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.416623) q[1];
sx q[1];
rz(-1.5048711) q[1];
sx q[1];
rz(-2.9030187) q[1];
x q[2];
rz(1.8252556) q[3];
sx q[3];
rz(-1.1512827) q[3];
sx q[3];
rz(2.6655243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7717379) q[2];
sx q[2];
rz(-2.3042046) q[2];
sx q[2];
rz(0.85956335) q[2];
rz(1.9178948) q[3];
sx q[3];
rz(-0.92646354) q[3];
sx q[3];
rz(0.74444509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0201465) q[0];
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
rz(1.0449833) q[3];
sx q[3];
rz(-0.079418728) q[3];
sx q[3];
rz(3.1374745) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
