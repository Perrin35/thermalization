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
rz(3.8745772) q[1];
sx q[1];
rz(12.193845) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46967888) q[0];
sx q[0];
rz(-2.1436474) q[0];
sx q[0];
rz(1.8125305) q[0];
rz(-pi) q[1];
rz(1.4650605) q[2];
sx q[2];
rz(-0.9423965) q[2];
sx q[2];
rz(-0.65199967) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8718611) q[1];
sx q[1];
rz(-1.4565804) q[1];
sx q[1];
rz(1.7608587) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7872693) q[3];
sx q[3];
rz(-1.6811922) q[3];
sx q[3];
rz(-2.20761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7131876) q[2];
sx q[2];
rz(-1.5822072) q[2];
sx q[2];
rz(2.2170846) q[2];
rz(-1.472578) q[3];
sx q[3];
rz(-1.2481097) q[3];
sx q[3];
rz(1.4991466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9532042) q[0];
sx q[0];
rz(-0.062347118) q[0];
sx q[0];
rz(-1.6054608) q[0];
rz(2.9470782) q[1];
sx q[1];
rz(-1.8201927) q[1];
sx q[1];
rz(0.054873437) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6560199) q[0];
sx q[0];
rz(-1.6632073) q[0];
sx q[0];
rz(-1.4569605) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61972159) q[2];
sx q[2];
rz(-1.1148858) q[2];
sx q[2];
rz(2.1330619) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3000852) q[1];
sx q[1];
rz(-1.9133139) q[1];
sx q[1];
rz(3.1373346) q[1];
rz(-pi) q[2];
rz(0.26344928) q[3];
sx q[3];
rz(-0.69034319) q[3];
sx q[3];
rz(3.1171947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7028971) q[2];
sx q[2];
rz(-0.56151152) q[2];
sx q[2];
rz(1.6595718) q[2];
rz(-2.1510018) q[3];
sx q[3];
rz(-1.4086658) q[3];
sx q[3];
rz(-1.3668485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7822587) q[0];
sx q[0];
rz(-0.5193091) q[0];
sx q[0];
rz(2.7666132) q[0];
rz(2.8938876) q[1];
sx q[1];
rz(-1.9539555) q[1];
sx q[1];
rz(1.3365655) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61993116) q[0];
sx q[0];
rz(-1.7921899) q[0];
sx q[0];
rz(3.0822166) q[0];
x q[1];
rz(-1.1281563) q[2];
sx q[2];
rz(-1.0269564) q[2];
sx q[2];
rz(1.3904872) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.8001576) q[1];
sx q[1];
rz(-0.86038024) q[1];
sx q[1];
rz(2.7374949) q[1];
x q[2];
rz(-0.7671719) q[3];
sx q[3];
rz(-1.8681521) q[3];
sx q[3];
rz(1.5945827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1027362) q[2];
sx q[2];
rz(-2.3114624) q[2];
sx q[2];
rz(2.1170199) q[2];
rz(1.864795) q[3];
sx q[3];
rz(-1.5494616) q[3];
sx q[3];
rz(-1.6259441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9643726) q[0];
sx q[0];
rz(-1.6765046) q[0];
sx q[0];
rz(-3.1080416) q[0];
rz(1.9112446) q[1];
sx q[1];
rz(-2.3760445) q[1];
sx q[1];
rz(-1.4039325) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0899635) q[0];
sx q[0];
rz(-1.7647867) q[0];
sx q[0];
rz(-0.35332638) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83547445) q[2];
sx q[2];
rz(-2.2611141) q[2];
sx q[2];
rz(0.26963216) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0701323) q[1];
sx q[1];
rz(-1.6825819) q[1];
sx q[1];
rz(1.2791355) q[1];
rz(-pi) q[2];
rz(-1.4675667) q[3];
sx q[3];
rz(-0.69403115) q[3];
sx q[3];
rz(-0.3895143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6491062) q[2];
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
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7810818) q[0];
sx q[0];
rz(-1.1163982) q[0];
sx q[0];
rz(1.0158585) q[0];
rz(-1.5083183) q[1];
sx q[1];
rz(-0.95817482) q[1];
sx q[1];
rz(3.1075081) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.236892) q[0];
sx q[0];
rz(-1.9634982) q[0];
sx q[0];
rz(1.7444201) q[0];
rz(-pi) q[1];
rz(1.836852) q[2];
sx q[2];
rz(-1.012946) q[2];
sx q[2];
rz(1.7709874) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8964256) q[1];
sx q[1];
rz(-2.8571269) q[1];
sx q[1];
rz(-1.5462272) q[1];
rz(-pi) q[2];
rz(2.199585) q[3];
sx q[3];
rz(-1.994547) q[3];
sx q[3];
rz(0.29867344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.74026996) q[2];
sx q[2];
rz(-0.57381845) q[2];
sx q[2];
rz(1.9704698) q[2];
rz(1.3327538) q[3];
sx q[3];
rz(-1.4274024) q[3];
sx q[3];
rz(3.0122421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.68833441) q[0];
sx q[0];
rz(-1.1880705) q[0];
sx q[0];
rz(-0.43584287) q[0];
rz(-1.1054989) q[1];
sx q[1];
rz(-1.2970129) q[1];
sx q[1];
rz(1.2058535) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77594513) q[0];
sx q[0];
rz(-1.8071113) q[0];
sx q[0];
rz(2.3274371) q[0];
rz(-0.016024307) q[2];
sx q[2];
rz(-1.4756225) q[2];
sx q[2];
rz(1.3118088) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7534415) q[1];
sx q[1];
rz(-0.51137629) q[1];
sx q[1];
rz(-0.053877342) q[1];
rz(-pi) q[2];
rz(-0.033127012) q[3];
sx q[3];
rz(-1.1713542) q[3];
sx q[3];
rz(-1.2101733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2445406) q[2];
sx q[2];
rz(-0.62770939) q[2];
sx q[2];
rz(-0.051606027) q[2];
rz(2.7339325) q[3];
sx q[3];
rz(-1.1838653) q[3];
sx q[3];
rz(-2.6962962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-0.62149858) q[0];
sx q[0];
rz(-0.61722732) q[0];
sx q[0];
rz(-0.67333418) q[0];
rz(-2.3576221) q[1];
sx q[1];
rz(-1.7242804) q[1];
sx q[1];
rz(0.53692445) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2818031) q[0];
sx q[0];
rz(-1.4872695) q[0];
sx q[0];
rz(3.0755755) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.827048) q[2];
sx q[2];
rz(-1.7125687) q[2];
sx q[2];
rz(1.8585376) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.463045) q[1];
sx q[1];
rz(-1.4765413) q[1];
sx q[1];
rz(2.3386392) q[1];
rz(2.5506006) q[3];
sx q[3];
rz(-1.6984852) q[3];
sx q[3];
rz(-1.5203116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9817104) q[2];
sx q[2];
rz(-1.8843702) q[2];
sx q[2];
rz(-2.9505777) q[2];
rz(0.28132176) q[3];
sx q[3];
rz(-1.1912991) q[3];
sx q[3];
rz(2.7991926) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0328338) q[0];
sx q[0];
rz(-2.7712951) q[0];
sx q[0];
rz(0.38761815) q[0];
rz(-0.11501137) q[1];
sx q[1];
rz(-1.4287881) q[1];
sx q[1];
rz(0.33755916) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96466366) q[0];
sx q[0];
rz(-1.7690072) q[0];
sx q[0];
rz(2.6508509) q[0];
x q[1];
rz(0.93162025) q[2];
sx q[2];
rz(-1.6281623) q[2];
sx q[2];
rz(-0.66040874) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2355289) q[1];
sx q[1];
rz(-1.7019094) q[1];
sx q[1];
rz(0.99793418) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0739325) q[3];
sx q[3];
rz(-2.3849871) q[3];
sx q[3];
rz(2.9697231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9541624) q[2];
sx q[2];
rz(-0.53415161) q[2];
sx q[2];
rz(-1.8743275) q[2];
rz(-2.3665442) q[3];
sx q[3];
rz(-1.5379484) q[3];
sx q[3];
rz(2.3601941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18701126) q[0];
sx q[0];
rz(-2.1033852) q[0];
sx q[0];
rz(-0.22098456) q[0];
rz(0.9221319) q[1];
sx q[1];
rz(-1.8716967) q[1];
sx q[1];
rz(1.0029213) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.636165) q[0];
sx q[0];
rz(-1.822585) q[0];
sx q[0];
rz(-1.4279143) q[0];
rz(-2.6384764) q[2];
sx q[2];
rz(-2.0965577) q[2];
sx q[2];
rz(2.2803277) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.50833118) q[1];
sx q[1];
rz(-1.08053) q[1];
sx q[1];
rz(2.0983178) q[1];
rz(-pi) q[2];
rz(-2.8691611) q[3];
sx q[3];
rz(-2.2386132) q[3];
sx q[3];
rz(1.9200793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6691436) q[2];
sx q[2];
rz(-0.74206918) q[2];
sx q[2];
rz(2.9830902) q[2];
rz(-2.6878099) q[3];
sx q[3];
rz(-2.3588389) q[3];
sx q[3];
rz(-2.2973072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.6145988) q[0];
sx q[0];
rz(-2.232593) q[0];
sx q[0];
rz(-0.19113834) q[0];
rz(0.29516164) q[1];
sx q[1];
rz(-2.252153) q[1];
sx q[1];
rz(-0.89231649) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3494209) q[0];
sx q[0];
rz(-0.76602174) q[0];
sx q[0];
rz(-1.3146521) q[0];
x q[1];
rz(0.093776137) q[2];
sx q[2];
rz(-0.9320335) q[2];
sx q[2];
rz(-0.057657777) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.72496966) q[1];
sx q[1];
rz(-1.6367216) q[1];
sx q[1];
rz(2.9030187) q[1];
rz(1.3163371) q[3];
sx q[3];
rz(-1.99031) q[3];
sx q[3];
rz(2.6655243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3698547) q[2];
sx q[2];
rz(-2.3042046) q[2];
sx q[2];
rz(0.85956335) q[2];
rz(-1.2236979) q[3];
sx q[3];
rz(-2.2151291) q[3];
sx q[3];
rz(-0.74444509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1214462) q[0];
sx q[0];
rz(-2.2059724) q[0];
sx q[0];
rz(-1.7511517) q[0];
rz(1.4795115) q[1];
sx q[1];
rz(-0.48304396) q[1];
sx q[1];
rz(1.2190291) q[1];
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