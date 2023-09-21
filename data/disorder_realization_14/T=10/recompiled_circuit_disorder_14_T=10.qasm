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
rz(-2.4560661315918) q[0];
sx q[0];
rz(5.89415469964082) q[0];
sx q[0];
rz(11.6827916860501) q[0];
rz(3.13186240196228) q[1];
sx q[1];
rz(4.59877303441102) q[1];
sx q[1];
rz(7.48143193720981) q[1];
cx q[1],q[0];
rz(1.04490125179291) q[0];
sx q[0];
rz(3.34310338099534) q[0];
sx q[0];
rz(10.4951212167661) q[0];
rz(-3.93688702583313) q[2];
sx q[2];
rz(4.47840598424012) q[2];
sx q[2];
rz(8.92249116896793) q[2];
cx q[2],q[1];
rz(3.42063570022583) q[1];
sx q[1];
rz(4.19940570195253) q[1];
sx q[1];
rz(5.90223667620822) q[1];
rz(1.52452301979065) q[3];
sx q[3];
rz(4.54178324540193) q[3];
sx q[3];
rz(7.89811310767337) q[3];
cx q[3],q[2];
rz(1.94642317295074) q[2];
sx q[2];
rz(5.30005064805085) q[2];
sx q[2];
rz(9.60612172483608) q[2];
rz(-0.261207342147827) q[3];
sx q[3];
rz(1.8758591731363) q[3];
sx q[3];
rz(10.1810881853025) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.83920347690582) q[0];
sx q[0];
rz(0.289724500971385) q[0];
sx q[0];
rz(12.95320556163) q[0];
rz(3.64398527145386) q[1];
sx q[1];
rz(4.11510440905625) q[1];
sx q[1];
rz(7.82505247592136) q[1];
cx q[1],q[0];
rz(-1.80973362922668) q[0];
sx q[0];
rz(3.88844952185685) q[0];
sx q[0];
rz(10.3684124708097) q[0];
rz(1.73506784439087) q[2];
sx q[2];
rz(4.26116088231141) q[2];
sx q[2];
rz(10.2363340616147) q[2];
cx q[2],q[1];
rz(-0.350404173135757) q[1];
sx q[1];
rz(4.70245459874208) q[1];
sx q[1];
rz(10.9598396778028) q[1];
rz(0.801170766353607) q[3];
sx q[3];
rz(5.65209475358064) q[3];
sx q[3];
rz(11.5593774080197) q[3];
cx q[3],q[2];
rz(-2.53884220123291) q[2];
sx q[2];
rz(4.0194241126352) q[2];
sx q[2];
rz(11.1517249107282) q[2];
rz(3.99192333221436) q[3];
sx q[3];
rz(2.70897048910195) q[3];
sx q[3];
rz(7.7414133310239) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.262811154127121) q[0];
sx q[0];
rz(4.67299905617768) q[0];
sx q[0];
rz(11.9561433553617) q[0];
rz(1.28948509693146) q[1];
sx q[1];
rz(4.1208423097902) q[1];
sx q[1];
rz(10.4167481422345) q[1];
cx q[1],q[0];
rz(1.00058555603027) q[0];
sx q[0];
rz(1.386454256373) q[0];
sx q[0];
rz(5.72113654612705) q[0];
rz(-0.902461647987366) q[2];
sx q[2];
rz(4.87195077736909) q[2];
sx q[2];
rz(10.147102034084) q[2];
cx q[2],q[1];
rz(-2.51326441764832) q[1];
sx q[1];
rz(-2.69184145132964) q[1];
sx q[1];
rz(8.06122086047336) q[1];
rz(-2.02815628051758) q[3];
sx q[3];
rz(4.4917507489496) q[3];
sx q[3];
rz(11.3070355415265) q[3];
cx q[3],q[2];
rz(-0.383018642663956) q[2];
sx q[2];
rz(6.44441691239411) q[2];
sx q[2];
rz(9.69306955336734) q[2];
rz(-0.394084364175797) q[3];
sx q[3];
rz(4.37252357800538) q[3];
sx q[3];
rz(9.23627064227268) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.85369336605072) q[0];
sx q[0];
rz(8.91108194191987) q[0];
sx q[0];
rz(9.82985652088329) q[0];
rz(-0.690080940723419) q[1];
sx q[1];
rz(4.29944792588288) q[1];
sx q[1];
rz(11.8685502767484) q[1];
cx q[1],q[0];
rz(-0.217012479901314) q[0];
sx q[0];
rz(6.37607184250886) q[0];
sx q[0];
rz(12.822272515289) q[0];
rz(-0.985445261001587) q[2];
sx q[2];
rz(5.65487638314302) q[2];
sx q[2];
rz(7.91060922145053) q[2];
cx q[2],q[1];
rz(-3.5039484500885) q[1];
sx q[1];
rz(3.66632875998551) q[1];
sx q[1];
rz(12.1929847955625) q[1];
rz(-0.60034054517746) q[3];
sx q[3];
rz(0.496490391092845) q[3];
sx q[3];
rz(10.553598022453) q[3];
cx q[3],q[2];
rz(2.4219388961792) q[2];
sx q[2];
rz(1.74393382866914) q[2];
sx q[2];
rz(7.91556558608218) q[2];
rz(0.404310673475266) q[3];
sx q[3];
rz(3.82410153945024) q[3];
sx q[3];
rz(11.0880944490354) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.88358354568481) q[0];
sx q[0];
rz(1.78593495686586) q[0];
sx q[0];
rz(11.5408324956815) q[0];
rz(0.571996629238129) q[1];
sx q[1];
rz(5.18881121476228) q[1];
sx q[1];
rz(16.3372835874478) q[1];
cx q[1],q[0];
rz(0.481101632118225) q[0];
sx q[0];
rz(4.41345504124696) q[0];
sx q[0];
rz(7.42403981684848) q[0];
rz(6.59724473953247) q[2];
sx q[2];
rz(0.700541885691234) q[2];
sx q[2];
rz(3.40899322032138) q[2];
cx q[2],q[1];
rz(0.242299005389214) q[1];
sx q[1];
rz(4.20119503338868) q[1];
sx q[1];
rz(6.15165874957248) q[1];
rz(1.30738067626953) q[3];
sx q[3];
rz(5.29728976090486) q[3];
sx q[3];
rz(11.0305084943692) q[3];
cx q[3],q[2];
rz(2.5489673614502) q[2];
sx q[2];
rz(5.91352716286714) q[2];
sx q[2];
rz(9.76712549328014) q[2];
rz(1.34583783149719) q[3];
sx q[3];
rz(1.44744089444215) q[3];
sx q[3];
rz(9.62079619466468) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.48900699615479) q[0];
sx q[0];
rz(5.04719558556611) q[0];
sx q[0];
rz(12.3811769246976) q[0];
rz(4.548095703125) q[1];
sx q[1];
rz(5.19226637681062) q[1];
sx q[1];
rz(8.0578036069791) q[1];
cx q[1],q[0];
rz(1.6706086397171) q[0];
sx q[0];
rz(2.14897224505479) q[0];
sx q[0];
rz(9.68956238626643) q[0];
rz(3.87158632278442) q[2];
sx q[2];
rz(4.94088652928407) q[2];
sx q[2];
rz(8.50153926610156) q[2];
cx q[2],q[1];
rz(-4.70252418518066) q[1];
sx q[1];
rz(-1.72261890570586) q[1];
sx q[1];
rz(7.94673857688113) q[1];
rz(-0.427545219659805) q[3];
sx q[3];
rz(2.09800902207429) q[3];
sx q[3];
rz(9.58118866979285) q[3];
cx q[3],q[2];
rz(2.8967444896698) q[2];
sx q[2];
rz(5.72977200348909) q[2];
sx q[2];
rz(10.3081539630811) q[2];
rz(1.7287437915802) q[3];
sx q[3];
rz(3.83404782612855) q[3];
sx q[3];
rz(10.2385651230733) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.890052318573) q[0];
sx q[0];
rz(3.244049692648) q[0];
sx q[0];
rz(10.7029936075132) q[0];
rz(-0.0378407686948776) q[1];
sx q[1];
rz(3.95692143042619) q[1];
sx q[1];
rz(10.8006583213727) q[1];
cx q[1],q[0];
rz(-1.02768874168396) q[0];
sx q[0];
rz(4.41922536690766) q[0];
sx q[0];
rz(11.8360595464627) q[0];
rz(4.65251922607422) q[2];
sx q[2];
rz(5.18522849877412) q[2];
sx q[2];
rz(14.1133737325589) q[2];
cx q[2],q[1];
rz(1.52328765392303) q[1];
sx q[1];
rz(4.59712079365785) q[1];
sx q[1];
rz(6.17699286936923) q[1];
rz(0.330261677503586) q[3];
sx q[3];
rz(1.62632790406282) q[3];
sx q[3];
rz(9.76430369018718) q[3];
cx q[3],q[2];
rz(-1.457479596138) q[2];
sx q[2];
rz(2.42396798928315) q[2];
sx q[2];
rz(8.69372240304157) q[2];
rz(0.11080064624548) q[3];
sx q[3];
rz(4.69807830651338) q[3];
sx q[3];
rz(13.2617470979612) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.596579968929291) q[0];
sx q[0];
rz(2.2667363007837) q[0];
sx q[0];
rz(10.1583422779958) q[0];
rz(-0.607975244522095) q[1];
sx q[1];
rz(1.94764509995515) q[1];
sx q[1];
rz(9.19048455952808) q[1];
cx q[1],q[0];
rz(3.74252820014954) q[0];
sx q[0];
rz(5.18747773964936) q[0];
sx q[0];
rz(10.7040588617246) q[0];
rz(7.12092161178589) q[2];
sx q[2];
rz(5.41524592240388) q[2];
sx q[2];
rz(7.91817162036105) q[2];
cx q[2],q[1];
rz(4.32582473754883) q[1];
sx q[1];
rz(0.824165256815501) q[1];
sx q[1];
rz(8.11312696932956) q[1];
rz(2.55572652816772) q[3];
sx q[3];
rz(1.28639355500276) q[3];
sx q[3];
rz(7.86667058467075) q[3];
cx q[3],q[2];
rz(0.126047879457474) q[2];
sx q[2];
rz(4.97390321095521) q[2];
sx q[2];
rz(10.8501212358396) q[2];
rz(-1.67832934856415) q[3];
sx q[3];
rz(8.64033094246919) q[3];
sx q[3];
rz(11.0707537889402) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.542068064212799) q[0];
sx q[0];
rz(0.33518377144868) q[0];
sx q[0];
rz(10.6181579589765) q[0];
rz(1.26061964035034) q[1];
sx q[1];
rz(4.50648644764955) q[1];
sx q[1];
rz(10.3695784568708) q[1];
cx q[1],q[0];
rz(-1.90766179561615) q[0];
sx q[0];
rz(2.23092702229554) q[0];
sx q[0];
rz(13.0614909887235) q[0];
rz(-0.630176901817322) q[2];
sx q[2];
rz(2.66532796819741) q[2];
sx q[2];
rz(8.47429052590534) q[2];
cx q[2],q[1];
rz(0.244308441877365) q[1];
sx q[1];
rz(5.95398131211335) q[1];
sx q[1];
rz(8.40795645713016) q[1];
rz(-0.355987995862961) q[3];
sx q[3];
rz(7.66892591317231) q[3];
sx q[3];
rz(10.627405858032) q[3];
cx q[3],q[2];
rz(-3.3580493927002) q[2];
sx q[2];
rz(5.21207323868806) q[2];
sx q[2];
rz(9.71298379301235) q[2];
rz(0.479734122753143) q[3];
sx q[3];
rz(5.23333683808381) q[3];
sx q[3];
rz(10.9327707052152) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.111868031322956) q[0];
sx q[0];
rz(6.55938163598115) q[0];
sx q[0];
rz(10.3377665638845) q[0];
rz(-0.374620139598846) q[1];
sx q[1];
rz(4.54500690300996) q[1];
sx q[1];
rz(10.3156896591107) q[1];
cx q[1],q[0];
rz(0.837153136730194) q[0];
sx q[0];
rz(7.60578313668305) q[0];
sx q[0];
rz(10.6499958991925) q[0];
rz(-0.393451809883118) q[2];
sx q[2];
rz(5.63273802598054) q[2];
sx q[2];
rz(9.75388393401309) q[2];
cx q[2],q[1];
rz(5.52156209945679) q[1];
sx q[1];
rz(2.74598106940324) q[1];
sx q[1];
rz(8.39068028926059) q[1];
rz(1.1924341917038) q[3];
sx q[3];
rz(3.34371955891187) q[3];
sx q[3];
rz(9.24128740131065) q[3];
cx q[3],q[2];
rz(3.37800335884094) q[2];
sx q[2];
rz(-0.858352510137014) q[2];
sx q[2];
rz(6.78479359149143) q[2];
rz(-1.28913986682892) q[3];
sx q[3];
rz(4.82983043988282) q[3];
sx q[3];
rz(8.97860623001262) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.298335313797) q[0];
sx q[0];
rz(2.18842175801332) q[0];
sx q[0];
rz(9.26600507496997) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(2.63801026344299) q[1];
sx q[1];
rz(1.50952378113801) q[1];
sx q[1];
rz(11.7505399942319) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(1.9737389087677) q[2];
sx q[2];
rz(1.38229945500428) q[2];
sx q[2];
rz(4.95237348078891) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-1.07907176017761) q[3];
sx q[3];
rz(3.79739579756791) q[3];
sx q[3];
rz(9.49401250331804) q[3];
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
