//Bring in dependencies
const express = require('express');
const exphbs = require('express-handlebars');
const path = require('path');
const mongoose = require('mongoose');
const bodyParser = require('body-parser');
const methodOverride = require('method-override');
const flash = require('connect-flash');
const session = require('express-session');
const passport = require('passport');
//Init app object
const app = express();

//MIDDLEWARE

//database config
//Check if using local dev or remote production database
const db = require('./config/database');

//mongoose - connect to db
mongoose.Promise = global.Promise;
mongoose.connect(db.mongoURL, {
  useMongoClient: true
})
  .then(() => console.log('MongoDB connected...'))
  .catch(err => console.log(err));

//handlebars
//This is just telling the system that we want to use the handlebars template engine
app.engine('handlebars', exphbs({
  //views dir will contain all our handlebar templates
  //The default layout will be set to main
  defaultLayout: 'main'
}));
app.set('view engine', 'handlebars');

//body parser middleware
app.use(bodyParser.urlencoded({encoded: false}));
app.use(bodyParser.json());

//method override middleware
app.use(methodOverride('_method'));

//session and flash middleware
app.use(session({
  secret: 'secret',
  resave: false,
  saveUninitialized: true
}));
app.use(flash());

//passport middleware
//This MUST be declared after the express session middleware
app.use(passport.initialize());
app.use(passport.session());

//GLOBAL VARS
app.use(function(request, response, next){
  //Flash messages
  response.locals.success_msg = request.flash('success_msg');
  response.locals.error_msg = request.flash('error_msg');
  response.locals.error = request.flash('error');
  //User info if logged in
  response.locals.user = request.user || null;
  next();
});

//SET STATIC DIRECTORY PATH
app.use(express.static(path.join(__dirname, 'public')));

//SETTING ROUTES
app.get('/', (request, response) => {
  const title = "Microbiology Data Portal"
  response.render('index', {
    title: title
  });
});

app.get('/vitek', (request, response) => {
  response.render('vitek');
});


//SET PORT
//You set a port for your listen method
const port = process.env.PORT || 5000;

//LISTEN METHOD
//Pass in the port value, and then pass in a callback function
app.listen(port, () => {
  console.log(`Server started on port ${port}`)
});
