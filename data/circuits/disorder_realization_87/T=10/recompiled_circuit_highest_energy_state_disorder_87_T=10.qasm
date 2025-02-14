OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.46621608734131) q[0];
sx q[0];
rz(5.33865466912324) q[0];
sx q[0];
rz(9.62617791294261) q[0];
rz(4.21387767791748) q[1];
sx q[1];
rz(3.90448513825471) q[1];
sx q[1];
rz(4.86210963725253) q[1];
cx q[1],q[0];
rz(1.38233768939972) q[0];
sx q[0];
rz(4.65755680401857) q[0];
sx q[0];
rz(11.9554140329282) q[0];
rz(2.46874094009399) q[2];
sx q[2];
rz(4.29320946534211) q[2];
sx q[2];
rz(10.7691111326139) q[2];
cx q[2],q[1];
rz(0.0464675426483154) q[1];
sx q[1];
rz(7.07058373292024) q[1];
sx q[1];
rz(7.2868919134061) q[1];
rz(-0.391864866018295) q[3];
sx q[3];
rz(-1.26149353186553) q[3];
sx q[3];
rz(5.88137838839694) q[3];
cx q[3],q[2];
rz(-0.0349155142903328) q[2];
sx q[2];
rz(2.33725288708741) q[2];
sx q[2];
rz(12.8454351186673) q[2];
rz(-1.25328242778778) q[3];
sx q[3];
rz(6.53294244607026) q[3];
sx q[3];
rz(12.4147703409116) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.61075532436371) q[0];
sx q[0];
rz(3.98881826003129) q[0];
sx q[0];
rz(9.6711647272031) q[0];
rz(1.94657444953918) q[1];
sx q[1];
rz(2.24768260319764) q[1];
sx q[1];
rz(11.0596986770551) q[1];
cx q[1],q[0];
rz(0.559132993221283) q[0];
sx q[0];
rz(4.46687856514985) q[0];
sx q[0];
rz(9.89041075705692) q[0];
rz(1.61368560791016) q[2];
sx q[2];
rz(4.40972331364686) q[2];
sx q[2];
rz(11.2654480695645) q[2];
cx q[2],q[1];
rz(1.82837021350861) q[1];
sx q[1];
rz(4.77015355427796) q[1];
sx q[1];
rz(9.78809458612605) q[1];
rz(3.66312313079834) q[3];
sx q[3];
rz(1.18809154828126) q[3];
sx q[3];
rz(8.97568229436084) q[3];
cx q[3],q[2];
rz(0.401114970445633) q[2];
sx q[2];
rz(4.11651453574235) q[2];
sx q[2];
rz(10.656766986839) q[2];
rz(2.039381980896) q[3];
sx q[3];
rz(3.14591111254925) q[3];
sx q[3];
rz(11.2298664808194) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.443228930234909) q[0];
sx q[0];
rz(2.44464275439317) q[0];
sx q[0];
rz(10.3326070070188) q[0];
rz(2.57461309432983) q[1];
sx q[1];
rz(7.26592508156831) q[1];
sx q[1];
rz(15.6051721334378) q[1];
cx q[1],q[0];
rz(2.58233213424683) q[0];
sx q[0];
rz(4.15304306347901) q[0];
sx q[0];
rz(12.8146710157315) q[0];
rz(3.24587917327881) q[2];
sx q[2];
rz(4.63181069691712) q[2];
sx q[2];
rz(13.8122215032498) q[2];
cx q[2],q[1];
rz(1.49900782108307) q[1];
sx q[1];
rz(8.68856206734712) q[1];
sx q[1];
rz(7.22903773783847) q[1];
rz(2.65850067138672) q[3];
sx q[3];
rz(2.53170409997041) q[3];
sx q[3];
rz(13.5816254377286) q[3];
cx q[3],q[2];
rz(3.42480707168579) q[2];
sx q[2];
rz(5.4745148738199) q[2];
sx q[2];
rz(14.0038184881131) q[2];
rz(0.392506748437881) q[3];
sx q[3];
rz(5.14812794526155) q[3];
sx q[3];
rz(9.12634745835468) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(4.23494672775269) q[0];
sx q[0];
rz(1.36273506482179) q[0];
sx q[0];
rz(10.6135732889096) q[0];
rz(-2.85547542572021) q[1];
sx q[1];
rz(5.66560259659822) q[1];
sx q[1];
rz(7.91245577334567) q[1];
cx q[1],q[0];
rz(0.855073690414429) q[0];
sx q[0];
rz(5.67980829079682) q[0];
sx q[0];
rz(12.062249159805) q[0];
rz(1.9173983335495) q[2];
sx q[2];
rz(0.798149736719676) q[2];
sx q[2];
rz(10.0236807823102) q[2];
cx q[2],q[1];
rz(-3.83922863006592) q[1];
sx q[1];
rz(5.65814915497834) q[1];
sx q[1];
rz(9.74548796414539) q[1];
rz(-1.78182888031006) q[3];
sx q[3];
rz(5.21012273629243) q[3];
sx q[3];
rz(5.71709296702548) q[3];
cx q[3],q[2];
rz(2.79012227058411) q[2];
sx q[2];
rz(4.22503057320649) q[2];
sx q[2];
rz(8.7477940082471) q[2];
rz(-1.24667680263519) q[3];
sx q[3];
rz(1.86919835408265) q[3];
sx q[3];
rz(11.0498912095945) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.93382728099823) q[0];
sx q[0];
rz(2.88498687942559) q[0];
sx q[0];
rz(9.44969142451092) q[0];
rz(-2.08615303039551) q[1];
sx q[1];
rz(5.29812303383882) q[1];
sx q[1];
rz(12.2829241514127) q[1];
cx q[1],q[0];
rz(-0.0585397779941559) q[0];
sx q[0];
rz(-1.31365665595) q[0];
sx q[0];
rz(12.8355688810269) q[0];
rz(-0.195853471755981) q[2];
sx q[2];
rz(1.69986596901948) q[2];
sx q[2];
rz(10.8013030052106) q[2];
cx q[2],q[1];
rz(3.53708934783936) q[1];
sx q[1];
rz(2.07141342957551) q[1];
sx q[1];
rz(7.81915352343723) q[1];
rz(-0.417342334985733) q[3];
sx q[3];
rz(3.85867634614045) q[3];
sx q[3];
rz(6.89101717471286) q[3];
cx q[3],q[2];
rz(2.82895660400391) q[2];
sx q[2];
rz(5.86572376092012) q[2];
sx q[2];
rz(15.3821301221769) q[2];
rz(-1.85879409313202) q[3];
sx q[3];
rz(1.15743401845033) q[3];
sx q[3];
rz(7.83078370093509) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.403898984193802) q[0];
sx q[0];
rz(1.65885451634461) q[0];
sx q[0];
rz(9.04981887935802) q[0];
rz(-2.66295790672302) q[1];
sx q[1];
rz(2.46919831831987) q[1];
sx q[1];
rz(12.9365188836972) q[1];
cx q[1],q[0];
rz(1.81068122386932) q[0];
sx q[0];
rz(6.25729385216767) q[0];
sx q[0];
rz(10.5276708364408) q[0];
rz(-2.95648145675659) q[2];
sx q[2];
rz(1.70776692231233) q[2];
sx q[2];
rz(10.4031299114148) q[2];
cx q[2],q[1];
rz(-0.00358414417132735) q[1];
sx q[1];
rz(4.04071262677247) q[1];
sx q[1];
rz(11.9629053831021) q[1];
rz(-0.00501359160989523) q[3];
sx q[3];
rz(5.14019146760041) q[3];
sx q[3];
rz(7.98930022715732) q[3];
cx q[3],q[2];
rz(-0.397735387086868) q[2];
sx q[2];
rz(3.33305944700772) q[2];
sx q[2];
rz(8.04371569155856) q[2];
rz(-1.04876518249512) q[3];
sx q[3];
rz(4.4865039904886) q[3];
sx q[3];
rz(11.9267587423246) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.253107070922852) q[0];
sx q[0];
rz(7.12237754662568) q[0];
sx q[0];
rz(11.6422147512357) q[0];
rz(0.916685163974762) q[1];
sx q[1];
rz(8.35856834252412) q[1];
sx q[1];
rz(7.68451223372623) q[1];
cx q[1],q[0];
rz(-1.22712540626526) q[0];
sx q[0];
rz(1.07733980019624) q[0];
sx q[0];
rz(9.49036097376748) q[0];
rz(-0.915650427341461) q[2];
sx q[2];
rz(7.71718040307099) q[2];
sx q[2];
rz(15.7667412519376) q[2];
cx q[2],q[1];
rz(-6.66255664825439) q[1];
sx q[1];
rz(4.55254689057405) q[1];
sx q[1];
rz(11.7367889642636) q[1];
rz(0.393518537282944) q[3];
sx q[3];
rz(5.7610434611612) q[3];
sx q[3];
rz(8.70920369624301) q[3];
cx q[3],q[2];
rz(-1.22831976413727) q[2];
sx q[2];
rz(4.82137301762635) q[2];
sx q[2];
rz(11.3413746118466) q[2];
rz(-0.00829588808119297) q[3];
sx q[3];
rz(6.27219239075715) q[3];
sx q[3];
rz(11.6992740392606) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.558157205581665) q[0];
sx q[0];
rz(5.45284429390962) q[0];
sx q[0];
rz(9.90504363774463) q[0];
rz(-3.28605389595032) q[1];
sx q[1];
rz(5.37552753289277) q[1];
sx q[1];
rz(5.41880748271152) q[1];
cx q[1],q[0];
rz(1.9821765422821) q[0];
sx q[0];
rz(4.6658457835489) q[0];
sx q[0];
rz(9.26103787719413) q[0];
rz(1.47612166404724) q[2];
sx q[2];
rz(1.60451582272584) q[2];
sx q[2];
rz(13.6911196470182) q[2];
cx q[2],q[1];
rz(2.49829530715942) q[1];
sx q[1];
rz(1.81132307847077) q[1];
sx q[1];
rz(9.97940067052051) q[1];
rz(0.279338926076889) q[3];
sx q[3];
rz(1.83508256276185) q[3];
sx q[3];
rz(9.29088813661739) q[3];
cx q[3],q[2];
rz(-3.3478410243988) q[2];
sx q[2];
rz(1.01030054886872) q[2];
sx q[2];
rz(6.91825888156101) q[2];
rz(-0.572799742221832) q[3];
sx q[3];
rz(4.49748578866059) q[3];
sx q[3];
rz(10.3001324891965) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.110058143734932) q[0];
sx q[0];
rz(2.36819711525971) q[0];
sx q[0];
rz(9.54904474168226) q[0];
rz(5.91793394088745) q[1];
sx q[1];
rz(4.5687839110666) q[1];
sx q[1];
rz(7.99323997496768) q[1];
cx q[1],q[0];
rz(1.30764031410217) q[0];
sx q[0];
rz(2.87882837851579) q[0];
sx q[0];
rz(7.10850689410373) q[0];
rz(-2.26341032981873) q[2];
sx q[2];
rz(8.64129033883149) q[2];
sx q[2];
rz(13.4543757200162) q[2];
cx q[2],q[1];
rz(4.00559043884277) q[1];
sx q[1];
rz(4.0528121312433) q[1];
sx q[1];
rz(12.9557916879575) q[1];
rz(-0.437677294015884) q[3];
sx q[3];
rz(5.19774761994416) q[3];
sx q[3];
rz(12.2066447496335) q[3];
cx q[3],q[2];
rz(-2.87538361549377) q[2];
sx q[2];
rz(5.36177721818025) q[2];
sx q[2];
rz(13.7391481161039) q[2];
rz(1.40699362754822) q[3];
sx q[3];
rz(4.9513717015558) q[3];
sx q[3];
rz(10.6394347906034) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.54919719696045) q[0];
sx q[0];
rz(4.27633622487123) q[0];
sx q[0];
rz(11.9717354535977) q[0];
rz(-1.00589621067047) q[1];
sx q[1];
rz(4.34400215943391) q[1];
sx q[1];
rz(13.4555611371915) q[1];
cx q[1],q[0];
rz(0.858146369457245) q[0];
sx q[0];
rz(3.36776825984056) q[0];
sx q[0];
rz(12.3385531663816) q[0];
rz(0.460365861654282) q[2];
sx q[2];
rz(1.6341007073694) q[2];
sx q[2];
rz(3.81640574931308) q[2];
cx q[2],q[1];
rz(0.883162498474121) q[1];
sx q[1];
rz(9.31296554406221) q[1];
sx q[1];
rz(8.8595602273862) q[1];
rz(-1.05440425872803) q[3];
sx q[3];
rz(4.74546566803987) q[3];
sx q[3];
rz(11.8113257646482) q[3];
cx q[3],q[2];
rz(0.951200067996979) q[2];
sx q[2];
rz(1.96596458752687) q[2];
sx q[2];
rz(6.82857367991611) q[2];
rz(2.17831802368164) q[3];
sx q[3];
rz(5.10721352894837) q[3];
sx q[3];
rz(7.74711809157535) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.26272630691528) q[0];
sx q[0];
rz(2.43383565743501) q[0];
sx q[0];
rz(8.15425417422458) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(1.44086909294128) q[1];
sx q[1];
rz(3.71725532610948) q[1];
sx q[1];
rz(10.4771106004636) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(2.03121113777161) q[2];
sx q[2];
rz(4.82748273213441) q[2];
sx q[2];
rz(9.60666067003413) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(0.0126286437734962) q[3];
sx q[3];
rz(3.83146235545213) q[3];
sx q[3];
rz(8.7034979224126) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
