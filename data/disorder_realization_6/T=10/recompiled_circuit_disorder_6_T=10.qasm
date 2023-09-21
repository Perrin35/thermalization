OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.525809943675995) q[0];
sx q[0];
rz(4.55941119988496) q[0];
sx q[0];
rz(8.8639080285947) q[0];
rz(4.25451183319092) q[1];
sx q[1];
rz(1.76340440114076) q[1];
sx q[1];
rz(7.49822852610751) q[1];
cx q[1],q[0];
rz(2.69988656044006) q[0];
sx q[0];
rz(2.92686483462388) q[0];
sx q[0];
rz(8.91570905446216) q[0];
rz(4.34808921813965) q[2];
sx q[2];
rz(3.83555004199082) q[2];
sx q[2];
rz(16.1348395109098) q[2];
cx q[2],q[1];
rz(2.34559679031372) q[1];
sx q[1];
rz(1.06182924111421) q[1];
sx q[1];
rz(8.23313972949191) q[1];
rz(4.56335353851318) q[3];
sx q[3];
rz(4.44970002968843) q[3];
sx q[3];
rz(10.9707933425824) q[3];
cx q[3],q[2];
rz(2.35401964187622) q[2];
sx q[2];
rz(5.33037963707978) q[2];
sx q[2];
rz(12.3832957506101) q[2];
rz(-0.377815276384354) q[3];
sx q[3];
rz(4.19038084347779) q[3];
sx q[3];
rz(5.98900029658481) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.297824442386627) q[0];
sx q[0];
rz(3.78631451924378) q[0];
sx q[0];
rz(6.20606801509067) q[0];
rz(3.48038959503174) q[1];
sx q[1];
rz(2.02705505688722) q[1];
sx q[1];
rz(7.8856478691022) q[1];
cx q[1],q[0];
rz(-1.21676433086395) q[0];
sx q[0];
rz(2.15688583453233) q[0];
sx q[0];
rz(10.9030778169553) q[0];
rz(1.33317673206329) q[2];
sx q[2];
rz(3.8944647033983) q[2];
sx q[2];
rz(10.6901416540067) q[2];
cx q[2],q[1];
rz(-0.512556076049805) q[1];
sx q[1];
rz(4.53814700444276) q[1];
sx q[1];
rz(8.93783879875346) q[1];
rz(0.235343962907791) q[3];
sx q[3];
rz(1.98415211041505) q[3];
sx q[3];
rz(11.2790544986646) q[3];
cx q[3],q[2];
rz(4.43765735626221) q[2];
sx q[2];
rz(5.00559225876863) q[2];
sx q[2];
rz(10.0832372665326) q[2];
rz(0.151308745145798) q[3];
sx q[3];
rz(4.16420439084107) q[3];
sx q[3];
rz(3.83650538920566) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.17005133628845) q[0];
sx q[0];
rz(3.92494604189927) q[0];
sx q[0];
rz(9.85788101553127) q[0];
rz(1.94948792457581) q[1];
sx q[1];
rz(1.21162179310853) q[1];
sx q[1];
rz(3.69694564341708) q[1];
cx q[1],q[0];
rz(0.992918908596039) q[0];
sx q[0];
rz(5.41613975365693) q[0];
sx q[0];
rz(10.0793010950009) q[0];
rz(5.0260534286499) q[2];
sx q[2];
rz(4.1694314797693) q[2];
sx q[2];
rz(13.6735296010892) q[2];
cx q[2],q[1];
rz(-2.80861163139343) q[1];
sx q[1];
rz(3.79738256533677) q[1];
sx q[1];
rz(12.8334927320401) q[1];
rz(-4.51102590560913) q[3];
sx q[3];
rz(7.31780496438081) q[3];
sx q[3];
rz(6.73823497294589) q[3];
cx q[3],q[2];
rz(-2.40425229072571) q[2];
sx q[2];
rz(8.63665405114228) q[2];
sx q[2];
rz(10.6753221511762) q[2];
rz(-0.244149699807167) q[3];
sx q[3];
rz(5.00096026261384) q[3];
sx q[3];
rz(7.97483465670749) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.260439485311508) q[0];
sx q[0];
rz(-0.457571593923024) q[0];
sx q[0];
rz(10.239579653732) q[0];
rz(-4.52097034454346) q[1];
sx q[1];
rz(3.49178627331788) q[1];
sx q[1];
rz(15.4527864217679) q[1];
cx q[1],q[0];
rz(-1.67386209964752) q[0];
sx q[0];
rz(2.6669589300924) q[0];
sx q[0];
rz(8.54466471671268) q[0];
rz(-1.13857126235962) q[2];
sx q[2];
rz(5.59722200234468) q[2];
sx q[2];
rz(11.4156596422116) q[2];
cx q[2],q[1];
rz(3.57082653045654) q[1];
sx q[1];
rz(5.9147893508249) q[1];
sx q[1];
rz(11.9476682901303) q[1];
rz(1.51977229118347) q[3];
sx q[3];
rz(2.09123662312562) q[3];
sx q[3];
rz(9.02114731668636) q[3];
cx q[3],q[2];
rz(4.03003025054932) q[2];
sx q[2];
rz(1.55068126519258) q[2];
sx q[2];
rz(9.59796223639652) q[2];
rz(-0.529824614524841) q[3];
sx q[3];
rz(6.13761487801606) q[3];
sx q[3];
rz(6.38551280497714) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.28168761730194) q[0];
sx q[0];
rz(4.79018596013124) q[0];
sx q[0];
rz(11.1904850959699) q[0];
rz(-1.86385226249695) q[1];
sx q[1];
rz(5.47100225289399) q[1];
sx q[1];
rz(9.368684398375) q[1];
cx q[1],q[0];
rz(1.84319055080414) q[0];
sx q[0];
rz(2.13923302491242) q[0];
sx q[0];
rz(12.1310765504758) q[0];
rz(1.223517537117) q[2];
sx q[2];
rz(4.27709999878938) q[2];
sx q[2];
rz(10.8842423915784) q[2];
cx q[2],q[1];
rz(1.54068148136139) q[1];
sx q[1];
rz(4.40964320500428) q[1];
sx q[1];
rz(9.40955144948467) q[1];
rz(-4.62314367294312) q[3];
sx q[3];
rz(2.12929943402345) q[3];
sx q[3];
rz(12.7974853277127) q[3];
cx q[3],q[2];
rz(0.000994801986962557) q[2];
sx q[2];
rz(5.36529413064057) q[2];
sx q[2];
rz(12.9985627889554) q[2];
rz(4.03578615188599) q[3];
sx q[3];
rz(4.24112811883027) q[3];
sx q[3];
rz(7.95865604876682) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.113192588090897) q[0];
sx q[0];
rz(2.25604614813859) q[0];
sx q[0];
rz(13.2139191389005) q[0];
rz(5.02127838134766) q[1];
sx q[1];
rz(4.81951919396455) q[1];
sx q[1];
rz(8.47028115986987) q[1];
cx q[1],q[0];
rz(-0.196767672896385) q[0];
sx q[0];
rz(4.55712810357148) q[0];
sx q[0];
rz(9.95814720391437) q[0];
rz(-1.60556674003601) q[2];
sx q[2];
rz(4.08209112484986) q[2];
sx q[2];
rz(6.73650262355014) q[2];
cx q[2],q[1];
rz(-0.336590141057968) q[1];
sx q[1];
rz(5.1878286917978) q[1];
sx q[1];
rz(11.5753686189572) q[1];
rz(0.291074186563492) q[3];
sx q[3];
rz(1.72725144227082) q[3];
sx q[3];
rz(8.48108092545673) q[3];
cx q[3],q[2];
rz(-0.592975616455078) q[2];
sx q[2];
rz(4.37506380875642) q[2];
sx q[2];
rz(8.38241419791385) q[2];
rz(-0.438676118850708) q[3];
sx q[3];
rz(4.19161930878694) q[3];
sx q[3];
rz(7.60126576422855) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.28386589884758) q[0];
sx q[0];
rz(6.05027952988679) q[0];
sx q[0];
rz(10.1679895281713) q[0];
rz(1.633953332901) q[1];
sx q[1];
rz(5.56329503853852) q[1];
sx q[1];
rz(8.81475559472247) q[1];
cx q[1],q[0];
rz(-1.2110720872879) q[0];
sx q[0];
rz(2.06063440640504) q[0];
sx q[0];
rz(9.13922542928859) q[0];
rz(3.79565763473511) q[2];
sx q[2];
rz(1.47872939904267) q[2];
sx q[2];
rz(9.66318452953502) q[2];
cx q[2],q[1];
rz(4.08337211608887) q[1];
sx q[1];
rz(3.80934390624101) q[1];
sx q[1];
rz(9.38432116284176) q[1];
rz(0.379294455051422) q[3];
sx q[3];
rz(4.52062788804109) q[3];
sx q[3];
rz(7.09667370318576) q[3];
cx q[3],q[2];
rz(0.336217761039734) q[2];
sx q[2];
rz(4.5839720090204) q[2];
sx q[2];
rz(10.3431812882344) q[2];
rz(-1.55047988891602) q[3];
sx q[3];
rz(2.1912626345926) q[3];
sx q[3];
rz(9.03587076663181) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.360887616872787) q[0];
sx q[0];
rz(3.81069943507249) q[0];
sx q[0];
rz(10.9382961749951) q[0];
rz(-0.529451191425323) q[1];
sx q[1];
rz(4.20831218560273) q[1];
sx q[1];
rz(11.8297853231351) q[1];
cx q[1],q[0];
rz(-0.373059570789337) q[0];
sx q[0];
rz(3.10964181100065) q[0];
sx q[0];
rz(9.77145159839793) q[0];
rz(1.46188008785248) q[2];
sx q[2];
rz(4.39185169537599) q[2];
sx q[2];
rz(9.46999526246592) q[2];
cx q[2],q[1];
rz(1.83839333057404) q[1];
sx q[1];
rz(2.03370776970918) q[1];
sx q[1];
rz(6.39708015917941) q[1];
rz(0.60968542098999) q[3];
sx q[3];
rz(0.529340656595775) q[3];
sx q[3];
rz(11.0562047719876) q[3];
cx q[3],q[2];
rz(4.17599773406982) q[2];
sx q[2];
rz(4.34907952149446) q[2];
sx q[2];
rz(5.59927341937228) q[2];
rz(1.22901391983032) q[3];
sx q[3];
rz(1.77142718632753) q[3];
sx q[3];
rz(10.8193182706754) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.19723153114319) q[0];
sx q[0];
rz(4.68434646924073) q[0];
sx q[0];
rz(9.23905349373027) q[0];
rz(-0.997052371501923) q[1];
sx q[1];
rz(4.40681448777253) q[1];
sx q[1];
rz(14.9631485700528) q[1];
cx q[1],q[0];
rz(-1.98415243625641) q[0];
sx q[0];
rz(2.30723485549028) q[0];
sx q[0];
rz(12.1249086618344) q[0];
rz(2.74061560630798) q[2];
sx q[2];
rz(2.32247606118257) q[2];
sx q[2];
rz(8.72372595071002) q[2];
cx q[2],q[1];
rz(-1.32651710510254) q[1];
sx q[1];
rz(2.72173667152459) q[1];
sx q[1];
rz(9.99876860379382) q[1];
rz(-0.392781615257263) q[3];
sx q[3];
rz(3.96894970734651) q[3];
sx q[3];
rz(7.10387823580905) q[3];
cx q[3],q[2];
rz(-1.03614377975464) q[2];
sx q[2];
rz(3.89732673962648) q[2];
sx q[2];
rz(11.3716896533887) q[2];
rz(0.996693253517151) q[3];
sx q[3];
rz(4.35765484173829) q[3];
sx q[3];
rz(11.5700125455777) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.88493347167969) q[0];
sx q[0];
rz(2.35345903237397) q[0];
sx q[0];
rz(9.82878101467296) q[0];
rz(-3.17271995544434) q[1];
sx q[1];
rz(4.6260017474466) q[1];
sx q[1];
rz(7.45413002967044) q[1];
cx q[1],q[0];
rz(0.0268965698778629) q[0];
sx q[0];
rz(5.85574284394319) q[0];
sx q[0];
rz(7.80170748233005) q[0];
rz(3.9367299079895) q[2];
sx q[2];
rz(4.51688245137269) q[2];
sx q[2];
rz(13.6385922193448) q[2];
cx q[2],q[1];
rz(3.5614161491394) q[1];
sx q[1];
rz(3.11739191797609) q[1];
sx q[1];
rz(6.6479267835538) q[1];
rz(1.14381766319275) q[3];
sx q[3];
rz(1.30377617676789) q[3];
sx q[3];
rz(9.51951314359113) q[3];
cx q[3],q[2];
rz(1.44604539871216) q[2];
sx q[2];
rz(4.92146936257417) q[2];
sx q[2];
rz(8.83281556367084) q[2];
rz(0.566360890865326) q[3];
sx q[3];
rz(3.30629593332345) q[3];
sx q[3];
rz(7.90094254016086) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.161501720547676) q[0];
sx q[0];
rz(0.704768093424388) q[0];
sx q[0];
rz(10.4593193292539) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(4.09388732910156) q[1];
sx q[1];
rz(5.69031730492646) q[1];
sx q[1];
rz(15.1047262906949) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-1.15522348880768) q[2];
sx q[2];
rz(2.45487496455247) q[2];
sx q[2];
rz(11.3926779985349) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-0.396210223436356) q[3];
sx q[3];
rz(5.57048335869844) q[3];
sx q[3];
rz(12.04216334819) q[3];
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