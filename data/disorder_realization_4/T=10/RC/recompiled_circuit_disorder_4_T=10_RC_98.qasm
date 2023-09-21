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
rz(0.92619196) q[0];
rz(1.9594833) q[1];
sx q[1];
rz(-0.73298454) q[1];
sx q[1];
rz(0.37252537) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46967888) q[0];
sx q[0];
rz(-2.1436474) q[0];
sx q[0];
rz(-1.3290622) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1442285) q[2];
sx q[2];
rz(-2.5055474) q[2];
sx q[2];
rz(0.47338212) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2791427) q[1];
sx q[1];
rz(-1.7596054) q[1];
sx q[1];
rz(3.0253009) q[1];
rz(-1.7872693) q[3];
sx q[3];
rz(-1.4604005) q[3];
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
rz(-0.92450809) q[2];
rz(-1.472578) q[3];
sx q[3];
rz(-1.893483) q[3];
sx q[3];
rz(-1.4991466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18838841) q[0];
sx q[0];
rz(-0.062347118) q[0];
sx q[0];
rz(1.6054608) q[0];
rz(2.9470782) q[1];
sx q[1];
rz(-1.3214) q[1];
sx q[1];
rz(3.0867192) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7354436) q[0];
sx q[0];
rz(-2.9950954) q[0];
sx q[0];
rz(2.2551401) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0286249) q[2];
sx q[2];
rz(-1.0222058) q[2];
sx q[2];
rz(2.2749535) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3000852) q[1];
sx q[1];
rz(-1.9133139) q[1];
sx q[1];
rz(3.1373346) q[1];
x q[2];
rz(0.67316405) q[3];
sx q[3];
rz(-1.4041956) q[3];
sx q[3];
rz(-1.8002321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.43869552) q[2];
sx q[2];
rz(-2.5800811) q[2];
sx q[2];
rz(-1.6595718) q[2];
rz(-2.1510018) q[3];
sx q[3];
rz(-1.7329268) q[3];
sx q[3];
rz(1.3668485) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7822587) q[0];
sx q[0];
rz(-0.5193091) q[0];
sx q[0];
rz(2.7666132) q[0];
rz(-2.8938876) q[1];
sx q[1];
rz(-1.1876371) q[1];
sx q[1];
rz(-1.8050271) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61993116) q[0];
sx q[0];
rz(-1.7921899) q[0];
sx q[0];
rz(-0.059376052) q[0];
x q[1];
rz(-1.1281563) q[2];
sx q[2];
rz(-2.1146362) q[2];
sx q[2];
rz(1.7511055) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6429236) q[1];
sx q[1];
rz(-1.8734697) q[1];
sx q[1];
rz(2.32294) q[1];
x q[2];
rz(-0.41575723) q[3];
sx q[3];
rz(-0.81167479) q[3];
sx q[3];
rz(2.8228001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0388564) q[2];
sx q[2];
rz(-2.3114624) q[2];
sx q[2];
rz(2.1170199) q[2];
rz(-1.864795) q[3];
sx q[3];
rz(-1.5494616) q[3];
sx q[3];
rz(-1.5156486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17722002) q[0];
sx q[0];
rz(-1.465088) q[0];
sx q[0];
rz(0.033551034) q[0];
rz(1.9112446) q[1];
sx q[1];
rz(-2.3760445) q[1];
sx q[1];
rz(-1.4039325) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0008154) q[0];
sx q[0];
rz(-0.40110943) q[0];
sx q[0];
rz(-2.6252069) q[0];
rz(-pi) q[1];
rz(-0.83547445) q[2];
sx q[2];
rz(-0.88047853) q[2];
sx q[2];
rz(-0.26963216) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6744068) q[1];
sx q[1];
rz(-1.8605839) q[1];
sx q[1];
rz(-3.0249216) q[1];
x q[2];
rz(0.87938829) q[3];
sx q[3];
rz(-1.6367568) q[3];
sx q[3];
rz(1.1018167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6491062) q[2];
sx q[2];
rz(-2.2968569) q[2];
sx q[2];
rz(2.9283294) q[2];
rz(0.02031859) q[3];
sx q[3];
rz(-1.3213986) q[3];
sx q[3];
rz(0.4030574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7810818) q[0];
sx q[0];
rz(-1.1163982) q[0];
sx q[0];
rz(1.0158585) q[0];
rz(1.6332743) q[1];
sx q[1];
rz(-2.1834178) q[1];
sx q[1];
rz(0.034084592) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.236892) q[0];
sx q[0];
rz(-1.1780945) q[0];
sx q[0];
rz(1.7444201) q[0];
rz(2.7427865) q[2];
sx q[2];
rz(-0.6119234) q[2];
sx q[2];
rz(-1.8460225) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.34921131) q[1];
sx q[1];
rz(-1.5776909) q[1];
sx q[1];
rz(1.2864119) q[1];
rz(-pi) q[2];
rz(0.50877737) q[3];
sx q[3];
rz(-2.1366589) q[3];
sx q[3];
rz(0.98154991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4013227) q[2];
sx q[2];
rz(-0.57381845) q[2];
sx q[2];
rz(-1.9704698) q[2];
rz(-1.3327538) q[3];
sx q[3];
rz(-1.7141902) q[3];
sx q[3];
rz(3.0122421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68833441) q[0];
sx q[0];
rz(-1.9535221) q[0];
sx q[0];
rz(2.7057498) q[0];
rz(-1.1054989) q[1];
sx q[1];
rz(-1.8445797) q[1];
sx q[1];
rz(1.9357392) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1036557) q[0];
sx q[0];
rz(-2.3559542) q[0];
sx q[0];
rz(1.2334137) q[0];
x q[1];
rz(1.4756104) q[2];
sx q[2];
rz(-1.5548445) q[2];
sx q[2];
rz(0.25746458) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0059433) q[1];
sx q[1];
rz(-1.597153) q[1];
sx q[1];
rz(2.6308358) q[1];
rz(-pi) q[2];
rz(1.6490963) q[3];
sx q[3];
rz(-2.7408528) q[3];
sx q[3];
rz(-2.0164255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.897052) q[2];
sx q[2];
rz(-0.62770939) q[2];
sx q[2];
rz(3.0899866) q[2];
rz(0.40766019) q[3];
sx q[3];
rz(-1.9577273) q[3];
sx q[3];
rz(0.44529644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62149858) q[0];
sx q[0];
rz(-2.5243653) q[0];
sx q[0];
rz(-0.67333418) q[0];
rz(0.78397059) q[1];
sx q[1];
rz(-1.7242804) q[1];
sx q[1];
rz(-2.6046682) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9519185) q[0];
sx q[0];
rz(-0.10641831) q[0];
sx q[0];
rz(-0.90344067) q[0];
x q[1];
rz(2.9951018) q[2];
sx q[2];
rz(-1.3171725) q[2];
sx q[2];
rz(-2.8908562) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1310136) q[1];
sx q[1];
rz(-2.3691633) q[1];
sx q[1];
rz(1.4355245) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5506006) q[3];
sx q[3];
rz(-1.6984852) q[3];
sx q[3];
rz(1.6212811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15988222) q[2];
sx q[2];
rz(-1.8843702) q[2];
sx q[2];
rz(-0.19101492) q[2];
rz(0.28132176) q[3];
sx q[3];
rz(-1.9502935) q[3];
sx q[3];
rz(0.34240001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1087588) q[0];
sx q[0];
rz(-2.7712951) q[0];
sx q[0];
rz(-2.7539745) q[0];
rz(0.11501137) q[1];
sx q[1];
rz(-1.4287881) q[1];
sx q[1];
rz(-0.33755916) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71096703) q[0];
sx q[0];
rz(-2.0511048) q[0];
sx q[0];
rz(-1.7947012) q[0];
x q[1];
rz(-0.07143306) q[2];
sx q[2];
rz(-2.2087503) q[2];
sx q[2];
rz(-0.86779867) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2355289) q[1];
sx q[1];
rz(-1.7019094) q[1];
sx q[1];
rz(-2.1436585) q[1];
rz(-1.5070595) q[3];
sx q[3];
rz(-0.81634854) q[3];
sx q[3];
rz(-0.07894978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.18743029) q[2];
sx q[2];
rz(-0.53415161) q[2];
sx q[2];
rz(1.2672651) q[2];
rz(0.77504843) q[3];
sx q[3];
rz(-1.6036443) q[3];
sx q[3];
rz(-2.3601941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
rz(2.2194608) q[1];
sx q[1];
rz(-1.8716967) q[1];
sx q[1];
rz(2.1386713) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1598338) q[0];
sx q[0];
rz(-2.8528385) q[0];
sx q[0];
rz(0.50555484) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2641364) q[2];
sx q[2];
rz(-0.71091953) q[2];
sx q[2];
rz(-0.029821776) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3468614) q[1];
sx q[1];
rz(-2.0309629) q[1];
sx q[1];
rz(-2.5882583) q[1];
x q[2];
rz(-0.27243154) q[3];
sx q[3];
rz(-0.90297943) q[3];
sx q[3];
rz(1.9200793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4724491) q[2];
sx q[2];
rz(-2.3995235) q[2];
sx q[2];
rz(-0.15850244) q[2];
rz(0.45378271) q[3];
sx q[3];
rz(-0.78275371) q[3];
sx q[3];
rz(2.2973072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6145988) q[0];
sx q[0];
rz(-0.90899962) q[0];
sx q[0];
rz(-2.9504543) q[0];
rz(0.29516164) q[1];
sx q[1];
rz(-0.8894397) q[1];
sx q[1];
rz(0.89231649) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96517262) q[0];
sx q[0];
rz(-1.3942379) q[0];
sx q[0];
rz(2.3202592) q[0];
rz(0.92992444) q[2];
sx q[2];
rz(-1.6460437) q[2];
sx q[2];
rz(1.6844695) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.82980624) q[1];
sx q[1];
rz(-1.3327507) q[1];
sx q[1];
rz(1.6386375) q[1];
rz(1.8252556) q[3];
sx q[3];
rz(-1.1512827) q[3];
sx q[3];
rz(2.6655243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3698547) q[2];
sx q[2];
rz(-2.3042046) q[2];
sx q[2];
rz(-0.85956335) q[2];
rz(1.9178948) q[3];
sx q[3];
rz(-0.92646354) q[3];
sx q[3];
rz(-2.3971476) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0201465) q[0];
sx q[0];
rz(-2.2059724) q[0];
sx q[0];
rz(-1.7511517) q[0];
rz(-1.6620811) q[1];
sx q[1];
rz(-0.48304396) q[1];
sx q[1];
rz(1.2190291) q[1];
rz(3.1142278) q[2];
sx q[2];
rz(-1.27925) q[2];
sx q[2];
rz(-1.7640511) q[2];
rz(0.039924351) q[3];
sx q[3];
rz(-1.502124) q[3];
sx q[3];
rz(0.52306642) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
