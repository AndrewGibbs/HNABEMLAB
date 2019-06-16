function [x,w] = encyclop_cub_square(d, N, v)
%function [x,w] = encyclop_cub_square(d, N, v)
%
%   Return a cubature formula of degree d and with N points. If there are
%   more possibilities, v=2,3,... will select them.
%   Source: Encyclopedia of Cubature Formulae (R. Cools)

% $Id$

if nargin == 2
    v = 1;
end

switch d
    case 3
        switch N
            case 4
                switch v
                    case 2
                        % ref [str71]
                        x1 = gen_fullysymmetric_aa(0.577350269189625764509148780501957);
                        w1 = gen_weights(x1, 1);
                        x = x1;
                        w = w1;
                    otherwise
                        % ref [str71]
                        x1 = gen_fullysymmetric_b0(0.816496580927726032732428024901963);
                        w1 = gen_weights(x1, 1);
                        x = x1;
                        w = w1;
                end
        end
    case 4
        switch N
            case 6
                % ref [wb86]
                x1 = gen_origin();
                w1 = gen_weights(x1, 1.14285714285714285714285714285714);
                x2 = gen_partialsymmetry_b0(0.966091783079295904914577610477984);
                w2 = gen_weights(x2, 0.439560439560439560439560439560439);
                x3 = gen_partialsymmetry(0.455603727836192842249584745989760, 0.851914653304600491301835640366669);
                w3 = gen_weights(x3, 0.566072207007532105264050459451812);
                x4 = gen_partialsymmetry(-0.731629951573134529368035491840613, 0.630912788976754026243413202993658);
                w4 = gen_weights(x4, 0.642719001783676685944740749339395);
                x = [x1; x2; x3; x4];
                w = [w1; w2; w3; w4];
        end
    case 5
        switch N
            case 7
                % ref [str71]
                x1 = gen_origin();
                w1 = gen_weights(x1, 1.14285714285714285714285714285714);
                x2 = gen_centralsymmetry(0.966091783079295904914577610477984, 0);
                w2 = gen_weights(x2, 0.317460317460317460317460317460317);
                x3 = gen_rectangularsymmetry(0.577350269189625764509148780501957, 0.774596669241483377035853079956479);
                w3 = gen_weights(x3, 0.555555555555555555555555555555555);
                x = [x1; x2; x3];
                w = [w1; w2; w3];
        end
    case 7
        switch N
            case 12
                switch v
                    case 2
                        % ref [hp77]
                        x1 = gen_rectangularsymmetry_b0(0.529422802042655325888559903752383);
                        w1 = gen_weights(x1, 0.635853883443279771824736563460042);
                        x2 = gen_rectangularsymmetry_a0(0.627041373780395317631777029507994);
                        w2 = gen_weights(x2, 0.590012715421030762966347426708615);
                        x3 = gen_rectangularsymmetry(0.917117822312770586255079396421675, 0.547931206828092323773641470296595);
                        w3 = gen_weights(x3, 0.213057211620949126506085869631973);
                        x4 = gen_rectangularsymmetry(0.611268766465328414395189713503039, 0.938843256658858304590191100766145);
                        w4 = gen_weights(x4, 0.174009488946895606098372135283697);
                        x = [x1; x2; x3; x4];
                        w = [w1; w2; w3; w4];
                    otherwise
                        % ref [str71]
                        x1 = gen_fullysymmetric_b0(0.925820099772551461566566776583999);
                        w1 = gen_weights(x1, 0.241975308641975308641975308641975);
                        x2 = gen_fullysymmetric_aa(0.380554433208315656379106359086394);
                        w2 = gen_weights(x2, 0.520592916667394457139919432046731);
                        x3 = gen_fullysymmetric_aa(0.805979782918598743707856181350744);
                        w3 = gen_weights(x3, 0.237431774690630234218105259311293);
                        x = [x1; x2; x3];
                        w = [w1; w2; w3];
                end
        end
    case 11
        switch N
            case 24
                x1 = gen_rotationalinvariant(0.982639223540855472952491497004009, 0.698076104549567564776469806174958);
                w1 = gen_weights(x1, 0.0480207633507238145627631759775806);
                x2 = gen_rotationalinvariant(0.825775835902963937302274585289940, 0.939486382816736907206432362169896);
                w2 = gen_weights(x2, 0.0660713291645505956736350808495464);
                x3 = gen_rotationalinvariant(0.188586138718641954600324568182745, 0.953539528201532015845004266823976);
                w3 = gen_weights(x3, 0.0973867773586681641961204397995472);
                x4 = gen_rotationalinvariant(0.812520548304813100489382581912299, 0.315623432915254195985609716402104);
                w4 = gen_weights(x4, 0.211736349998948600503931661356261);
                x5 = gen_rotationalinvariant(0.525320250364547762341631887140024, 0.712001913075336306549065895123759);
                w5 = gen_weights(x5, 0.225626061728863387403158016208490);
                x6 = gen_rotationalinvariant(0.0416580719120223682735468045377018, 0.424847248848669250615430111511957);
                w6 = gen_weights(x6, 0.351158718398245437660391625808574);
                x = [x1; x2; x3; x4; x5; x6];
                w = [w1; w2; w3; w4; w5; w6];
        end
    case 17
        switch N
            case 60
                x1 = gen_fullysymmetric_b0(0.989353074512600491226435837616855);
                w1 = gen_weights(x1, 0.0206149159199909598849346111352593);
                x2 = gen_fullysymmetric_b0(0.376285207157973294407548611720249);
                w2 = gen_weights(x2, 0.128025716179909834718790578266771);
                x3 = gen_fullysymmetric_aa(0.978848279262233116070071435993900);
                w3 = gen_weights(x3, 5.51173953403189052015467719393391e-3);
                x4 = gen_fullysymmetric_aa(0.885794729164116128906454358065004);
                w4 = gen_weights(x4, 0.0392077124571418804955998450600411);
                x5 = gen_fullysymmetric_aa(0.171756123838348174694383891493585);
                w5 = gen_weights(x5, 0.0763969450798633024005048543017064);
                x6 = gen_fullysymmetric(0.590499273806002413351516521008016, 0.319505036634573945676789474303272);
                w6 = gen_weights(x6, 0.141513729949972459244609749412573);
                x7 = gen_fullysymmetric(0.799079131916863255810133334567146, 0.597972451929457380262801755996745);
                w7 = gen_weights(x7, 0.0839032793637976020823276661920304);
                x8 = gen_fullysymmetric(0.803743962958744711785449825494367, 0.0583444817765505296843531563292007);
                w8 = gen_weights(x8, 0.0603941636496845460309035890917432);
                x9 = gen_fullysymmetric(0.936506276127494781742088348171630, 0.347386316166202674785997791701337);
                w9 = gen_weights(x9, 0.0573877529692126951214158298337486);
                x10 = gen_fullysymmetric(0.981321179805452294975910440430711, 0.706000287798646119181926220623163);
                w10 = gen_weights(x10, 0.0219225594818637635107508824910476);
                x = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10];
                w = [w1; w2; w3; w4; w5; w6; w7; w8; w9; w10];
        end
end


function z = gen_weights(x1, y)
z = ones(size(x1,1),1) * y;

function z = gen_origin()
z = [0 0];

function z = gen_partialsymmetry_b0(a)
b = 0;
z = [a b];

function z = gen_partialsymmetry(a, b)
z = [a b; a -b];

function z = gen_fullysymmetric_b0(a)
b = 0;
z = [a b; -a b; b a; b -a];

function z = gen_fullysymmetric_a0(b)
a = 0;
z = [a b; a -b; b a; -b a];

function z = gen_fullysymmetric_aa(a)
z = [a a; a -a; -a a; -a -a];

function z = gen_fullysymmetric(a, b)
z = [a b; -a b; a -b; -a -b; b a; -b a; b -a; -b -a];

function z = gen_centralsymmetry(a, b)
z = [a b; -a -b];

function z = gen_rectangularsymmetry(a, b)
z = [a b; -a b; a -b; -a -b];

function z = gen_rectangularsymmetry_a0(b)
a = 0;
z = [a b; a -b];

function z = gen_rectangularsymmetry_b0(a)
b = 0;
z = [a b; -a b];

function z = gen_rotationalinvariant(a, b)
z = [a b; -a -b; -b a; b -a];

