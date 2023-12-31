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
rz(1.9594833) q[1];
sx q[1];
rz(-0.73298454) q[1];
sx q[1];
rz(-2.7690673) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9076219) q[0];
sx q[0];
rz(-1.7733493) q[0];
sx q[0];
rz(0.58637182) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5105272) q[2];
sx q[2];
rz(-1.4853146) q[2];
sx q[2];
rz(-2.1604872) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8358742) q[1];
sx q[1];
rz(-0.22138518) q[1];
sx q[1];
rz(2.1165044) q[1];
rz(-pi) q[2];
rz(1.7872693) q[3];
sx q[3];
rz(-1.6811922) q[3];
sx q[3];
rz(0.9339827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.42840502) q[2];
sx q[2];
rz(-1.5822072) q[2];
sx q[2];
rz(-2.2170846) q[2];
rz(1.6690147) q[3];
sx q[3];
rz(-1.2481097) q[3];
sx q[3];
rz(-1.6424461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9532042) q[0];
sx q[0];
rz(-3.0792455) q[0];
sx q[0];
rz(-1.6054608) q[0];
rz(-0.19451441) q[1];
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
rz(-0.48557278) q[0];
sx q[0];
rz(-1.6632073) q[0];
sx q[0];
rz(-1.4569605) q[0];
x q[1];
rz(1.0286249) q[2];
sx q[2];
rz(-1.0222058) q[2];
sx q[2];
rz(-0.86663914) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.84150746) q[1];
sx q[1];
rz(-1.2282787) q[1];
sx q[1];
rz(3.1373346) q[1];
rz(-1.7826471) q[3];
sx q[3];
rz(-0.90862521) q[3];
sx q[3];
rz(-0.3609095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.43869552) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7822587) q[0];
sx q[0];
rz(-2.6222836) q[0];
sx q[0];
rz(-2.7666132) q[0];
rz(2.8938876) q[1];
sx q[1];
rz(-1.9539555) q[1];
sx q[1];
rz(1.3365655) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61993116) q[0];
sx q[0];
rz(-1.3494028) q[0];
sx q[0];
rz(3.0822166) q[0];
rz(-0.61632421) q[2];
sx q[2];
rz(-0.68683544) q[2];
sx q[2];
rz(0.64885215) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3414351) q[1];
sx q[1];
rz(-0.86038024) q[1];
sx q[1];
rz(0.40409778) q[1];
rz(-2.7258354) q[3];
sx q[3];
rz(-0.81167479) q[3];
sx q[3];
rz(0.31879253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0388564) q[2];
sx q[2];
rz(-0.83013022) q[2];
sx q[2];
rz(-1.0245727) q[2];
rz(-1.2767977) q[3];
sx q[3];
rz(-1.5494616) q[3];
sx q[3];
rz(1.5156486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9643726) q[0];
sx q[0];
rz(-1.465088) q[0];
sx q[0];
rz(3.1080416) q[0];
rz(-1.230348) q[1];
sx q[1];
rz(-0.76554811) q[1];
sx q[1];
rz(1.4039325) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0008154) q[0];
sx q[0];
rz(-0.40110943) q[0];
sx q[0];
rz(0.51638575) q[0];
x q[1];
rz(0.68217512) q[2];
sx q[2];
rz(-0.96207679) q[2];
sx q[2];
rz(-1.9145554) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.4671859) q[1];
sx q[1];
rz(-1.2810088) q[1];
sx q[1];
rz(-3.0249216) q[1];
rz(-3.0560533) q[3];
sx q[3];
rz(-0.88118689) q[3];
sx q[3];
rz(-2.6181108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4924865) q[2];
sx q[2];
rz(-2.2968569) q[2];
sx q[2];
rz(0.21326324) q[2];
rz(-0.02031859) q[3];
sx q[3];
rz(-1.3213986) q[3];
sx q[3];
rz(2.7385353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.236892) q[0];
sx q[0];
rz(-1.9634982) q[0];
sx q[0];
rz(1.3971726) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5675315) q[2];
sx q[2];
rz(-1.3458369) q[2];
sx q[2];
rz(0.056919295) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2451671) q[1];
sx q[1];
rz(-0.28446576) q[1];
sx q[1];
rz(1.5953654) q[1];
rz(-pi) q[2];
rz(2.6328153) q[3];
sx q[3];
rz(-2.1366589) q[3];
sx q[3];
rz(2.1600427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.74026996) q[2];
sx q[2];
rz(-2.5677742) q[2];
sx q[2];
rz(1.1711228) q[2];
rz(-1.8088388) q[3];
sx q[3];
rz(-1.4274024) q[3];
sx q[3];
rz(3.0122421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4532582) q[0];
sx q[0];
rz(-1.9535221) q[0];
sx q[0];
rz(-2.7057498) q[0];
rz(-2.0360937) q[1];
sx q[1];
rz(-1.2970129) q[1];
sx q[1];
rz(1.9357392) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1036557) q[0];
sx q[0];
rz(-2.3559542) q[0];
sx q[0];
rz(-1.2334137) q[0];
rz(3.1255683) q[2];
sx q[2];
rz(-1.4756225) q[2];
sx q[2];
rz(1.3118088) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1356493) q[1];
sx q[1];
rz(-1.5444396) q[1];
sx q[1];
rz(0.51075682) q[1];
x q[2];
rz(3.1084656) q[3];
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
rz(1.897052) q[2];
sx q[2];
rz(-0.62770939) q[2];
sx q[2];
rz(0.051606027) q[2];
rz(2.7339325) q[3];
sx q[3];
rz(-1.1838653) q[3];
sx q[3];
rz(0.44529644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5200941) q[0];
sx q[0];
rz(-2.5243653) q[0];
sx q[0];
rz(-2.4682585) q[0];
rz(2.3576221) q[1];
sx q[1];
rz(-1.4173123) q[1];
sx q[1];
rz(0.53692445) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1896742) q[0];
sx q[0];
rz(-0.10641831) q[0];
sx q[0];
rz(-2.238152) q[0];
x q[1];
rz(-0.14649086) q[2];
sx q[2];
rz(-1.8244201) q[2];
sx q[2];
rz(-0.25073642) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9432224) q[1];
sx q[1];
rz(-2.3343625) q[1];
sx q[1];
rz(-0.13065773) q[1];
rz(-pi) q[2];
rz(-1.4173996) q[3];
sx q[3];
rz(-2.1563357) q[3];
sx q[3];
rz(-0.034754001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.15988222) q[2];
sx q[2];
rz(-1.2572224) q[2];
sx q[2];
rz(0.19101492) q[2];
rz(0.28132176) q[3];
sx q[3];
rz(-1.1912991) q[3];
sx q[3];
rz(2.7991926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1087588) q[0];
sx q[0];
rz(-0.37029752) q[0];
sx q[0];
rz(-0.38761815) q[0];
rz(-0.11501137) q[1];
sx q[1];
rz(-1.7128046) q[1];
sx q[1];
rz(-0.33755916) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25305504) q[0];
sx q[0];
rz(-0.5262143) q[0];
sx q[0];
rz(-2.7387268) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.07143306) q[2];
sx q[2];
rz(-0.93284235) q[2];
sx q[2];
rz(-2.273794) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.90606373) q[1];
sx q[1];
rz(-1.4396832) q[1];
sx q[1];
rz(0.99793418) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5070595) q[3];
sx q[3];
rz(-2.3252441) q[3];
sx q[3];
rz(0.07894978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9541624) q[2];
sx q[2];
rz(-2.607441) q[2];
sx q[2];
rz(1.2672651) q[2];
rz(2.3665442) q[3];
sx q[3];
rz(-1.5379484) q[3];
sx q[3];
rz(0.78139853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.9545814) q[0];
sx q[0];
rz(-2.1033852) q[0];
sx q[0];
rz(-2.9206081) q[0];
rz(-2.2194608) q[1];
sx q[1];
rz(-1.8716967) q[1];
sx q[1];
rz(-2.1386713) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-2.2641364) q[2];
sx q[2];
rz(-2.4306731) q[2];
sx q[2];
rz(3.1117709) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.79473125) q[1];
sx q[1];
rz(-1.1106297) q[1];
sx q[1];
rz(-2.5882583) q[1];
rz(-pi) q[2];
rz(0.88463155) q[3];
sx q[3];
rz(-1.3579218) q[3];
sx q[3];
rz(-0.5205982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4724491) q[2];
sx q[2];
rz(-2.3995235) q[2];
sx q[2];
rz(-0.15850244) q[2];
rz(2.6878099) q[3];
sx q[3];
rz(-0.78275371) q[3];
sx q[3];
rz(0.84428549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
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
rz(-0.8894397) q[1];
sx q[1];
rz(-0.89231649) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6979881) q[0];
sx q[0];
rz(-2.3058878) q[0];
sx q[0];
rz(2.9025335) q[0];
x q[1];
rz(3.0478165) q[2];
sx q[2];
rz(-0.9320335) q[2];
sx q[2];
rz(-3.0839349) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1103507) q[1];
sx q[1];
rz(-0.24734766) q[1];
sx q[1];
rz(-0.27242839) q[1];
x q[2];
rz(0.51384135) q[3];
sx q[3];
rz(-2.6548879) q[3];
sx q[3];
rz(-1.0443618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3698547) q[2];
sx q[2];
rz(-2.3042046) q[2];
sx q[2];
rz(0.85956335) q[2];
rz(1.9178948) q[3];
sx q[3];
rz(-0.92646354) q[3];
sx q[3];
rz(-2.3971476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
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
rz(1.6620811) q[1];
sx q[1];
rz(-2.6585487) q[1];
sx q[1];
rz(-1.9225635) q[1];
rz(-1.8624458) q[2];
sx q[2];
rz(-1.5970061) q[2];
sx q[2];
rz(2.9562052) q[2];
rz(1.5020694) q[3];
sx q[3];
rz(-1.5309661) q[3];
sx q[3];
rz(-1.0504709) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
