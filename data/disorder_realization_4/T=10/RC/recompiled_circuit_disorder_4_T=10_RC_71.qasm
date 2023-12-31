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
rz(0.37252537) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042779048) q[0];
sx q[0];
rz(-2.5251237) q[0];
sx q[0];
rz(0.35538506) q[0];
rz(-2.5105272) q[2];
sx q[2];
rz(-1.6562781) q[2];
sx q[2];
rz(-0.98110547) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.26973154) q[1];
sx q[1];
rz(-1.4565804) q[1];
sx q[1];
rz(-1.3807339) q[1];
rz(-pi) q[2];
rz(2.0472237) q[3];
sx q[3];
rz(-0.24260394) q[3];
sx q[3];
rz(-2.0403595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42840502) q[2];
sx q[2];
rz(-1.5593854) q[2];
sx q[2];
rz(-0.92450809) q[2];
rz(-1.6690147) q[3];
sx q[3];
rz(-1.893483) q[3];
sx q[3];
rz(1.4991466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.8201927) q[1];
sx q[1];
rz(-3.0867192) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0957735) q[0];
sx q[0];
rz(-1.4574483) q[0];
sx q[0];
rz(-0.093009526) q[0];
x q[1];
rz(-0.70116455) q[2];
sx q[2];
rz(-0.75116457) q[2];
sx q[2];
rz(-3.1322111) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.82882985) q[1];
sx q[1];
rz(-0.34254303) q[1];
sx q[1];
rz(1.5588552) q[1];
rz(-pi) q[2];
rz(-0.26344928) q[3];
sx q[3];
rz(-0.69034319) q[3];
sx q[3];
rz(0.024397959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7028971) q[2];
sx q[2];
rz(-2.5800811) q[2];
sx q[2];
rz(1.6595718) q[2];
rz(-0.99059087) q[3];
sx q[3];
rz(-1.4086658) q[3];
sx q[3];
rz(-1.7747442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3593339) q[0];
sx q[0];
rz(-2.6222836) q[0];
sx q[0];
rz(2.7666132) q[0];
rz(2.8938876) q[1];
sx q[1];
rz(-1.9539555) q[1];
sx q[1];
rz(-1.8050271) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5216615) q[0];
sx q[0];
rz(-1.7921899) q[0];
sx q[0];
rz(-0.059376052) q[0];
rz(-pi) q[1];
rz(2.0134363) q[2];
sx q[2];
rz(-1.0269564) q[2];
sx q[2];
rz(-1.7511055) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3805483) q[1];
sx q[1];
rz(-2.3420463) q[1];
sx q[1];
rz(1.9995081) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9732479) q[3];
sx q[3];
rz(-2.2964722) q[3];
sx q[3];
rz(-2.890051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1027362) q[2];
sx q[2];
rz(-2.3114624) q[2];
sx q[2];
rz(1.0245727) q[2];
rz(1.2767977) q[3];
sx q[3];
rz(-1.5921311) q[3];
sx q[3];
rz(1.5156486) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9643726) q[0];
sx q[0];
rz(-1.465088) q[0];
sx q[0];
rz(-0.033551034) q[0];
rz(1.230348) q[1];
sx q[1];
rz(-0.76554811) q[1];
sx q[1];
rz(1.7376602) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0516292) q[0];
sx q[0];
rz(-1.3768059) q[0];
sx q[0];
rz(-0.35332638) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83547445) q[2];
sx q[2];
rz(-2.2611141) q[2];
sx q[2];
rz(-2.8719605) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.4671859) q[1];
sx q[1];
rz(-1.8605839) q[1];
sx q[1];
rz(-3.0249216) q[1];
x q[2];
rz(2.2622044) q[3];
sx q[3];
rz(-1.5048358) q[3];
sx q[3];
rz(1.1018167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6491062) q[2];
sx q[2];
rz(-0.84473574) q[2];
sx q[2];
rz(2.9283294) q[2];
rz(0.02031859) q[3];
sx q[3];
rz(-1.820194) q[3];
sx q[3];
rz(-0.4030574) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7810818) q[0];
sx q[0];
rz(-2.0251944) q[0];
sx q[0];
rz(2.1257341) q[0];
rz(1.5083183) q[1];
sx q[1];
rz(-0.95817482) q[1];
sx q[1];
rz(0.034084592) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9047006) q[0];
sx q[0];
rz(-1.9634982) q[0];
sx q[0];
rz(-1.3971726) q[0];
x q[1];
rz(-0.39880619) q[2];
sx q[2];
rz(-0.6119234) q[2];
sx q[2];
rz(-1.8460225) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7923813) q[1];
sx q[1];
rz(-1.5776909) q[1];
sx q[1];
rz(-1.2864119) q[1];
x q[2];
rz(-2.2250416) q[3];
sx q[3];
rz(-0.74186462) q[3];
sx q[3];
rz(1.3548917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68833441) q[0];
sx q[0];
rz(-1.9535221) q[0];
sx q[0];
rz(0.43584287) q[0];
rz(-1.1054989) q[1];
sx q[1];
rz(-1.2970129) q[1];
sx q[1];
rz(-1.9357392) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77594513) q[0];
sx q[0];
rz(-1.3344814) q[0];
sx q[0];
rz(-0.81415557) q[0];
rz(-pi) q[1];
x q[1];
rz(0.016024307) q[2];
sx q[2];
rz(-1.6659701) q[2];
sx q[2];
rz(-1.8297838) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0059433) q[1];
sx q[1];
rz(-1.5444396) q[1];
sx q[1];
rz(0.51075682) q[1];
rz(1.1711575) q[3];
sx q[3];
rz(-1.6013147) q[3];
sx q[3];
rz(-2.7680824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2445406) q[2];
sx q[2];
rz(-2.5138833) q[2];
sx q[2];
rz(0.051606027) q[2];
rz(-0.40766019) q[3];
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
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62149858) q[0];
sx q[0];
rz(-2.5243653) q[0];
sx q[0];
rz(0.67333418) q[0];
rz(-2.3576221) q[1];
sx q[1];
rz(-1.7242804) q[1];
sx q[1];
rz(-2.6046682) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1896742) q[0];
sx q[0];
rz(-3.0351743) q[0];
sx q[0];
rz(-0.90344067) q[0];
rz(-pi) q[1];
x q[1];
rz(0.14649086) q[2];
sx q[2];
rz(-1.8244201) q[2];
sx q[2];
rz(-2.8908562) q[2];
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
rz(2.9817104) q[2];
sx q[2];
rz(-1.2572224) q[2];
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
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0328338) q[0];
sx q[0];
rz(-2.7712951) q[0];
sx q[0];
rz(0.38761815) q[0];
rz(0.11501137) q[1];
sx q[1];
rz(-1.4287881) q[1];
sx q[1];
rz(2.8040335) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.176929) q[0];
sx q[0];
rz(-1.7690072) q[0];
sx q[0];
rz(-2.6508509) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2099724) q[2];
sx q[2];
rz(-1.5134303) q[2];
sx q[2];
rz(0.66040874) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2355289) q[1];
sx q[1];
rz(-1.4396832) q[1];
sx q[1];
rz(2.1436585) q[1];
rz(-pi) q[2];
rz(3.0739325) q[3];
sx q[3];
rz(-0.75660556) q[3];
sx q[3];
rz(2.9697231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.18743029) q[2];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-1.2698959) q[1];
sx q[1];
rz(-1.0029213) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10119469) q[0];
sx q[0];
rz(-1.4324491) q[0];
sx q[0];
rz(-2.8873216) q[0];
rz(-pi) q[1];
rz(2.1557751) q[2];
sx q[2];
rz(-1.1406116) q[2];
sx q[2];
rz(-2.1625724) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.79473125) q[1];
sx q[1];
rz(-2.0309629) q[1];
sx q[1];
rz(-0.55333432) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88463155) q[3];
sx q[3];
rz(-1.3579218) q[3];
sx q[3];
rz(-0.5205982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4724491) q[2];
sx q[2];
rz(-2.3995235) q[2];
sx q[2];
rz(-2.9830902) q[2];
rz(2.6878099) q[3];
sx q[3];
rz(-2.3588389) q[3];
sx q[3];
rz(2.2973072) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6145988) q[0];
sx q[0];
rz(-2.232593) q[0];
sx q[0];
rz(-0.19113834) q[0];
rz(-0.29516164) q[1];
sx q[1];
rz(-2.252153) q[1];
sx q[1];
rz(0.89231649) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79217171) q[0];
sx q[0];
rz(-0.76602174) q[0];
sx q[0];
rz(1.3146521) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4453663) q[2];
sx q[2];
rz(-0.64465603) q[2];
sx q[2];
rz(0.21412011) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3117864) q[1];
sx q[1];
rz(-1.808842) q[1];
sx q[1];
rz(-1.5029552) q[1];
rz(-1.3163371) q[3];
sx q[3];
rz(-1.1512827) q[3];
sx q[3];
rz(-0.47606836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7717379) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1214462) q[0];
sx q[0];
rz(-2.2059724) q[0];
sx q[0];
rz(-1.7511517) q[0];
rz(-1.4795115) q[1];
sx q[1];
rz(-2.6585487) q[1];
sx q[1];
rz(-1.9225635) q[1];
rz(-1.66172) q[2];
sx q[2];
rz(-2.8488013) q[2];
sx q[2];
rz(-1.6691096) q[2];
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
